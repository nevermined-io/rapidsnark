#ifndef ZKEY_UTILS_PLONK_H
#define ZKEY_UTILS_PLONK_H

#include <gmp.h>
#include <memory>

#include "binfile_utils.hpp"

namespace ZKeyUtilsPlonk {

    class Header {


    public:
        u_int32_t n8q;
        mpz_t qPrime;
        u_int32_t n8r;
        mpz_t rPrime;

        u_int32_t nVars;
        u_int32_t nPublic;
        u_int32_t domainSize;
        u_int64_t nAdditions;
        u_int64_t nConstrains;

        mpz_t k1;
        mpz_t k2;

        void *Qm;
        void *Ql;
        void *Qr;
        void *Qo;
        void *Qc;
        void *S1;
        void *S2;
        void *S3;
        void *X_2;

        Header();
        ~Header();
    };

    Header *loadHeader(BinFileUtils::BinFile *f);
}

#endif // ZKEY_UTILS_H
