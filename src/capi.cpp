
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
#include "utils.hpp"
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

typedef void (*ItFunc)(Circom_CalcWit *ctx, int idx, json val);

void iterateArr(Circom_CalcWit *ctx, int o, Circom_Sizes sizes, json jarr, ItFunc f) {
  if (!jarr.is_array()) {
    assert((sizes[0] == 1)&&(sizes[1] == 0));
    f(ctx, o, jarr);
  } else {
    int n = sizes[0] / sizes[1];
    for (int i=0; i<n; i++) {
      iterateArr(ctx, o + i*sizes[1], sizes+1, jarr[i], f);
    }
  }
}

void itFunc(Circom_CalcWit *ctx, int o, json val) {

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
        handle_error("Invalid JSON type");
    }

    Fr_str2element (&v, s.c_str());

    ctx->setSignal(0, 0, o, &v);
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
      int o;
      try {
        o = ctx->getSignalOffset(0, h);
      } catch (std::runtime_error e) {
        std::ostringstream errStrStream;
        errStrStream << "Error loadin variable: " << it.key() << "\n" << e.what();
        throw std::runtime_error(errStrStream.str() );
      }
      Circom_Sizes sizes = ctx->getSignalSizes(0, h);
      iterateArr(ctx, o, sizes, it.value(), itFunc);
    }
}

#define ADJ_P(a) *((void **)&a) = (void *)(((char *)circuit)+ (uint64_t)(a))

Circom_Circuit *loadCircuit(std::string const &datFileName) {
    Circom_Circuit *circuitF;
    Circom_Circuit *circuit;

    int fd;
    struct stat sb;

    fd = open(datFileName.c_str(), O_RDONLY);
    if (fd == -1) {
        std::cout << ".dat file not found: " << datFileName << "\n";
        throw std::system_error(errno, std::generic_category(), "open");
    }

    if (fstat(fd, &sb) == -1) {         /* To obtain file size */
        throw std::system_error(errno, std::generic_category(), "fstat");
    }

    circuitF = (Circom_Circuit *)mmap(NULL, sb.st_size, PROT_READ , MAP_PRIVATE, fd, 0);
    close(fd);

    circuit = (Circom_Circuit *)malloc(sb.st_size);
    memcpy((void *)circuit, (void *)circuitF, sb.st_size);

    munmap(circuitF, sb.st_size);

    ADJ_P(circuit->wit2sig);
    ADJ_P(circuit->components);
    ADJ_P(circuit->mapIsInput);
    ADJ_P(circuit->constants);
    ADJ_P(circuit->P);
    ADJ_P(circuit->componentEntries);

    for (int i=0; i<circuit->NComponents; i++) {
        ADJ_P(circuit->components[i].hashTable);
        ADJ_P(circuit->components[i].entries);
        circuit->components[i].fn = _functionTable[  (uint64_t)circuit->components[i].fn];
    }

    for (int i=0; i<circuit->NComponentEntries; i++) {
        ADJ_P(circuit->componentEntries[i].sizes);
    }

    return circuit;
}

void writeOutBin(Circom_CalcWit *ctx, std::string filename) {
    FILE *write_ptr;

    write_ptr = fopen(filename.c_str(),"wb");

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

    u32 nVars = circuit->NVars;
    fwrite(&nVars, 4, 1, write_ptr);


    // Data
    u32 idSection2 = 2;
    fwrite(&idSection2, 4, 1, write_ptr);

    u64 idSection2length = (u64)n8*(u64)circuit->NVars;
    fwrite(&idSection2length, 8, 1, write_ptr);

    FrElement v;

    for (int i=0;i<circuit->NVars;i++) {
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
        std::cerr << wtnsFilename << " " << inFilename << "\n";

        /*
        dta->calcwit->reset();
        loadJson(dta->calcwit, inFilename);
        dta->calcwit->join();
        writeOutBin(dta->calcwit, wtnsFilename);
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
        */
       return 0;
    }

/*
    int main() {
        void *ptr = make("keytransfer.zkey", "keytransfer.dat");
        fullprove(ptr, "tmp.wtns", "input.json");
    }
*/
}

