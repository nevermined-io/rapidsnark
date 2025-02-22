
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <gmp.h>
#include <memory>
#include <stdexcept>
#include <nlohmann/json.hpp>

#include <alt_bn128.hpp>
#include "binfile_utils.hpp"
#include "zkey_utils_plonk.hpp"
#include "wtns_utils.hpp"
#include "plonk.hpp"

#include "calcwit.hpp"
#include "circom.hpp"
#include <sstream>
#include <string>
#include <iomanip>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>

using json = nlohmann::json;

Circom_Circuit *circuit;

#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)
struct my_data {
    Plonk::Prover<AltBn128::Engine> *prover;
    Circom_CalcWit *calcwit;
};


void json2FrElements (json val, std::vector<FrElement> & vval){
  if (!val.is_array()) {
    FrElement v;
    std::string s;
    if (val.is_string()) {
        s = val.get<std::string>();
    } else if (val.is_number()) {
        double vd = val.get<double>();
        std::stringstream stream;
        stream << std::fixed << std::setprecision(0) << vd;
        s = stream.str();
    } else {
        throw new std::runtime_error("Invalid JSON type");
    }
    Fr_str2element (&v, s.c_str());
    vval.push_back(v);
  } else {
    for (uint i = 0; i < val.size(); i++) {
      json2FrElements (val[i], vval);
    }
  }
}


void loadJson(Circom_CalcWit *ctx, std::string filename) {
  std::ifstream inStream(filename);
  json j;
  inStream >> j;
  
  u64 nItems = j.size();
  // printf("Items : %llu\n",nItems);
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    // std::cout << it.key() << " => " << it.value() << '\n';
    u64 h = fnv1a(it.key());
    std::vector<FrElement> v;
    json2FrElements(it.value(),v);
    for (uint i = 0; i<v.size(); i++){
      try {
        // std::cout << it.key() << "," << i << " => " << Fr_element2str(&(v[i])) << '\n';
        ctx->setInputSignal(h,i,v[i]);
      } catch (std::runtime_error e) {
        std::ostringstream errStrStream;
        errStrStream << "Error loading variable: " << it.key() << "\n" << e.what();
        throw std::runtime_error(errStrStream.str() );
      }
    }
  }
}

Circom_Circuit* loadCircuit(std::string const &datFileName) {
    Circom_Circuit *circuit = new Circom_Circuit;

    int fd;
    struct stat sb;

    fd = open(datFileName.c_str(), O_RDONLY);
    if (fd == -1) {
        std::cout << ".dat file not found: " << datFileName << "\n";
        throw std::system_error(errno, std::generic_category(), "open");
    }
    
    if (fstat(fd, &sb) == -1) {          /* To obtain file size */
        throw std::system_error(errno, std::generic_category(), "fstat");
    }

    u8* bdata = (u8*)mmap(NULL, sb.st_size, PROT_READ , MAP_PRIVATE, fd, 0);
    close(fd);

    circuit->InputHashMap = new HashSignalInfo[get_size_of_input_hashmap()];
    uint dsize = get_size_of_input_hashmap()*sizeof(HashSignalInfo);
    memcpy((void *)(circuit->InputHashMap), (void *)bdata, dsize);

    circuit->witness2SignalList = new u64[get_size_of_witness()];
    uint inisize = dsize;    
    dsize = get_size_of_witness()*sizeof(u64);
    memcpy((void *)(circuit->witness2SignalList), (void *)(bdata+inisize), dsize);

    circuit->circuitConstants = new FrElement[get_size_of_constants()];
    if (get_size_of_constants()>0) {
      inisize += dsize;
      dsize = get_size_of_constants()*sizeof(FrElement);
      memcpy((void *)(circuit->circuitConstants), (void *)(bdata+inisize), dsize);
    }

    std::map<u32,IODefPair> templateInsId2IOSignalInfo1;
    if (get_size_of_io_map()>0) {
      u32 index[get_size_of_io_map()];
      inisize += dsize;
      dsize = get_size_of_io_map()*sizeof(u32);
      memcpy((void *)index, (void *)(bdata+inisize), dsize);
      inisize += dsize;
      assert(inisize % sizeof(u32) == 0);    
      assert(sb.st_size % sizeof(u32) == 0);
      u32 dataiomap[(sb.st_size-inisize)/sizeof(u32)];
      memcpy((void *)dataiomap, (void *)(bdata+inisize), sb.st_size-inisize);
      u32* pu32 = dataiomap;

      for (int i = 0; i < get_size_of_io_map(); i++) {
        u32 n = *pu32;
        IODefPair p;
        p.len = n;
        IODef defs[n];
        pu32 += 1;
        for (u32 j = 0; j <n; j++){
          defs[j].offset=*pu32;
          u32 len = *(pu32+1);
          defs[j].len = len;
          defs[j].lengths = new u32[len];
          memcpy((void *)defs[j].lengths,(void *)(pu32+2),len*sizeof(u32));
          pu32 += len + 2;
        }
        p.defs = (IODef*)calloc(10, sizeof(IODef));
        for (u32 j = 0; j < p.len; j++){
          p.defs[j] = defs[j];
        }
        templateInsId2IOSignalInfo1[index[i]] = p;
      }
    }
    circuit->templateInsId2IOSignalInfo = move(templateInsId2IOSignalInfo1);
    
    munmap(bdata, sb.st_size);
    
    return circuit;
}

void writeBinWitness(Circom_CalcWit *ctx, std::string wtnsFileName) {
    FILE *write_ptr;

    write_ptr = fopen(wtnsFileName.c_str(),"wb");

    fwrite("wtns", 4, 1, write_ptr);

    u32 version = 2;
    fwrite(&version, 4, 1, write_ptr);

    u32 nSections = 2;
    fwrite(&nSections, 4, 1, write_ptr);

    // Header
    u32 idSection1 = 1;
    fwrite(&idSection1, 4, 1, write_ptr);

    u32 n8 = Fr_N64*8;

    u64 idSection1length = 8 + n8;
    fwrite(&idSection1length, 8, 1, write_ptr);

    fwrite(&n8, 4, 1, write_ptr);

    fwrite(Fr_q.longVal, Fr_N64*8, 1, write_ptr);

    uint Nwtns = get_size_of_witness();
    
    u32 nVars = (u32)Nwtns;
    fwrite(&nVars, 4, 1, write_ptr);

    // Data
    u32 idSection2 = 2;
    fwrite(&idSection2, 4, 1, write_ptr);
    
    u64 idSection2length = (u64)n8*(u64)Nwtns;
    fwrite(&idSection2length, 8, 1, write_ptr);

    FrElement v;

    for (int i=0;i<Nwtns;i++) {
        ctx->getWitness(i, &v);
        Fr_toLongNormal(&v, &v);
        fwrite(v.longVal, Fr_N64*8, 1, write_ptr);
    }
    fclose(write_ptr);
}

extern "C" {
    void *make(const char *_zkey, const char *_dat) {
        Logger::getInstance()->enableConsoleLogging();
        Logger::getInstance()->updateLogLevel(DISABLE_LOG);
        // Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);
        std::string zkeyFilename = _zkey;
        std::string datFilename = _dat;

        circuit = loadCircuit(datFilename);

        // open output
        Circom_CalcWit *ctx = new Circom_CalcWit(circuit);

        auto zkey = BinFileUtils::openExisting(zkeyFilename, "zkey", 1);
        auto zkeyHeader = ZKeyUtilsPlonk::loadHeader(zkey);

        auto prover = Plonk::makeProver<AltBn128::Engine>(
            zkeyHeader->nVars,
            zkeyHeader->nPublic,
            zkeyHeader->domainSize,
            zkeyHeader->nAdditions,
            zkeyHeader->nConstrains,
            zkeyHeader->Qm,
            zkeyHeader->Ql,
            zkeyHeader->Qr,
            zkeyHeader->Qo,
            zkeyHeader->Qc,
            zkeyHeader->S1,
            zkeyHeader->S2,
            zkeyHeader->S3,
            zkeyHeader->X_2,
            zkeyHeader->n8r,
            zkeyHeader->n8q,
            zkey->getSectionData(3),    // Additions
            zkey->getSectionData(4),    // Amap
            zkey->getSectionData(5),    // Bmap
            zkey->getSectionData(6),    // Cmap
            zkey->getSectionData(14),   // Ptau
            zkeyHeader->k1,
            zkeyHeader->k2,
            zkey->getSectionData(12),   // Sigma
            zkey->getSectionData(7),   // Qm
            zkey->getSectionData(8),   // Ql
            zkey->getSectionData(9),   // Qr
            zkey->getSectionData(10),   // Qo
            zkey->getSectionData(11),   // Qc
            zkey->getSectionData(13)    // polData
        );
        my_data *dta = new my_data;
        dta->prover = prover;
        // std::cerr << "Made prover " << (uint64_t)dta << "\n";
        dta->calcwit = ctx;
        return (void*)dta;
    }

    char *fullprove(void *ptr, const char *_wtns, const char *_in) {
        my_data *dta = (my_data*)ptr;
        std::string wtnsFilename = _wtns;
        std::string inFilename = _in;

        // TODO: make sure reset and join aren't needed anymore
        dta->calcwit->reset();
        loadJson(dta->calcwit, inFilename);
        // dta->calcwit->join();
        writeBinWitness(dta->calcwit, wtnsFilename);
        // std::cerr << "wrote binary\n";

        auto wtns = BinFileUtils::openExisting(wtnsFilename, "wtns", 2);
        auto wtnsHeader = WtnsUtils::loadHeader(wtns);

        AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);
        // std::cerr << "Going to prove " << (uint64_t)ptr << "\n";
        auto proof = dta->prover->prove(wtnsData);
        // std::cerr << "Proof: " << proof << "\n";
        char *result = new char[proof.size()+1];
        for (int i = 0; i < proof.size(); i++) {
            result[i] = proof[i];
        }
        result[proof.size()] = 0;
        return result;
    }

/*
    int main() {
        void *ptr = make("keytransfer.zkey", "keytransfer.dat");
        fullprove(ptr, "tmp.wtns", "input.json");
    }
*/
}

