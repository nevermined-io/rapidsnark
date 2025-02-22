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

using json = nlohmann::json;

#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)

int main(int argc, char **argv) {

    if (argc != 5) {
        std::cerr << "Invalid number of parameters:\n";
        std::cerr << "Usage: prove <circuit.zkey> <witness.wtns> <proof.json> <public.json>\n";
        return -1;
    }

    mpz_t altBbn128r;

    mpz_init(altBbn128r);
    mpz_set_str(altBbn128r, "21888242871839275222246405745257275088548364400416034343698204186575808495617", 10);

    Logger::getInstance()->enableConsoleLogging();
    Logger::getInstance()->updateLogLevel(LOG_LEVEL_DEBUG);

    try {
        std::string zkeyFilename = argv[1];
        std::string wtnsFilename = argv[2];
        std::string proofFilename = argv[3];
        std::string publicFilename = argv[4];

        auto zkey = BinFileUtils::openExisting(zkeyFilename, "zkey", 1);
        auto zkeyHeader = ZKeyUtilsPlonk::loadHeader(zkey);

        std::string proofStr;
        if (mpz_cmp(zkeyHeader->rPrime, altBbn128r) != 0) {
            throw std::invalid_argument( "zkey curve not supported" );
        }

        auto wtns = BinFileUtils::openExisting(wtnsFilename, "wtns", 2);
        auto wtnsHeader = WtnsUtils::loadHeader(wtns);

        if (mpz_cmp(wtnsHeader->prime, altBbn128r) != 0) {
            throw std::invalid_argument( "different wtns curve" );
        }

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
        AltBn128::FrElement *wtnsData = (AltBn128::FrElement *)wtns->getSectionData(2);
        auto proof = prover->prove(wtnsData);

        std::ofstream proofFile;
        proofFile.open (proofFilename);
        proofFile << proof;
        proofFile.close();

        std::cerr << "Proof: " << proof << "\n";

        std::ofstream publicFile;
        publicFile.open (publicFilename);

        json jsonPublic;
        AltBn128::FrElement aux;
        for (int i=1; i<=zkeyHeader->nPublic; i++) {
            AltBn128::Fr.toMontgomery(aux, wtnsData[i]);
            jsonPublic.push_back(AltBn128::Fr.toString(aux));
        }

        publicFile << jsonPublic;
        publicFile.close();

    } catch (std::exception& e) {
        mpz_clear(altBbn128r);
        std::cerr << e.what() << '\n';
        return -1;
    }

    mpz_clear(altBbn128r);
    exit(EXIT_SUCCESS);
}
