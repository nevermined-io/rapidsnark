#include <stdexcept>

#include "zkey_utils_plonk.hpp"

namespace ZKeyUtilsPlonk {


Header::Header() {
}

Header::~Header() {
    mpz_clear(qPrime);
    mpz_clear(rPrime);
}

Header *loadHeader(BinFileUtils::BinFile *f) {
    auto h = new Header();

    f->startReadSection(1);
    uint32_t protocol = f->readU32LE();
    if (protocol != 2) {
        throw std::invalid_argument( "zkey file is not plonk" );
    }
    f->endReadSection();

    f->startReadSection(2);

    h->n8q = f->readU32LE();
    mpz_init(h->qPrime);
    mpz_import(h->qPrime, h->n8q, -1, 1, -1, 0, f->read(h->n8q));

    h->n8r = f->readU32LE();
    mpz_init(h->rPrime);
    mpz_import(h->rPrime, h->n8r , -1, 1, -1, 0, f->read(h->n8r));

    h->nVars = f->readU32LE();
    h->nPublic = f->readU32LE();
    h->domainSize = f->readU32LE();
    h->nAdditions = f->readU32LE();
    h->nConstrains = f->readU32LE();

    mpz_init(h->k1);
    mpz_import(h->k1, h->n8r , -1, 1, -1, 0, f->read(h->n8r));

    mpz_init(h->k2);
    mpz_import(h->k2, h->n8r , -1, 1, -1, 0, f->read(h->n8r));

    h->Qm = f->read(h->n8q*2);
    h->Ql = f->read(h->n8q*2);
    h->Qr = f->read(h->n8q*2);
    h->Qo = f->read(h->n8q*2);
    h->Qc = f->read(h->n8q*2);
    h->S1 = f->read(h->n8q*2);
    h->S2 = f->read(h->n8q*2);
    h->S3 = f->read(h->n8q*2);
    h->X_2 = f->read(h->n8q*4);
    f->endReadSection();

    return h;
}

} // namespace

