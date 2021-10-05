#include <sodium.h>

#include "logger.hpp"

using namespace CPlusPlusLogging;

namespace Plonk {


u_int32_t readUInt32(void *buf, uint32_t pos) {
    return *((u_int32_t *)((u_int64_t)buf + pos));
}


template <typename Engine>
std::unique_ptr<Proof<Engine>> Prover<Engine>::prove(typename Engine::FrElement *_wtns) {
    wtns = _wtns;
    LOG_DEBUG(E.f2.toString(X_2.x).c_str());
    wtns[0] = E.fr.zero();

    // Calculate additions
    // typename Engine::FrElement *internalWitness = new typename Engine::FrElement[n8r*nAdditions];
    uint32_t sSum = 8+n8r*2;

    typename Engine::FrElement r;
    typename Engine::FrElement aux1;
    typename Engine::FrElement aux2;
    for (uint32_t i = 0; i < nAdditions; i++) {
        uint32_t ai = readUInt32(additions, i*sSum);
        uint32_t bi = readUInt32(additions, i*sSum+4);
        // LOG_DEBUG(std::to_string(ai));
        // LOG_DEBUG(std::to_string(bi));
        typename Engine::FrElement& ac = *(typename Engine::FrElement *)(additions+i*sSum+8);
        typename Engine::FrElement& bc = *(typename Engine::FrElement *)(additions+i*sSum+8+n8r);
        typename Engine::FrElement& aw = getWitness(ai);
        typename Engine::FrElement& bw = getWitness(bi);
        // LOG_DEBUG("got witnesses");
        // LOG_DEBUG(E.fr.toString(aw).c_str());
        // LOG_DEBUG(E.fr.toString(bw).c_str());

        E.fr.mul(aux1, ac, aw);
        // LOG_DEBUG(E.fr.toString(aux1).c_str());
        // LOG_DEBUG("mul");
        E.fr.mul(aux2, bc, bw);
        // LOG_DEBUG("add");
        E.fr.add(r, aux1, aux2);
        internalWitness[i] = r;
    }
    LOG_DEBUG("Computed additions");
    
    // Round1

    // Build ABC
    typename Engine::FrElement *A = new typename Engine::FrElement[domainSize];
    typename Engine::FrElement *B = new typename Engine::FrElement[domainSize];
    typename Engine::FrElement *C = new typename Engine::FrElement[domainSize];

    for (uint32_t i = 0; i < domainSize; i++) {
        A[i] = E.fr.zero();
    }

    for (uint32_t i = 0; i < nConstrain; i++) {
        uint32_t iA = readUInt32(Amap, i*4);
        A[i] = getWitness(iA);
        if (i == 14) {
            LOG_DEBUG("Check witness 14");
            LOG_DEBUG(std::to_string(iA));
        }
        uint32_t iB = readUInt32(Bmap, i*4);
        B[i] = getWitness(iB);
        uint32_t iC = readUInt32(Cmap, i*4);
        C[i] = getWitness(iC);
    }

    LOG_DEBUG(E.fr.toString(A[0]).c_str());
    LOG_DEBUG("Check bar 14");
    LOG_DEBUG(E.fr.toString(A[14]).c_str());

    for (uint32_t i = 0; i < nConstrain; i++) {
        E.fr.toMontgomery(A[i], A[i]);
        E.fr.toMontgomery(B[i], B[i]);
        E.fr.toMontgomery(C[i], C[i]);
    }

    LOG_DEBUG(E.fr.toString(A[0]).c_str());

    // Random elements
    typename Engine::FrElement *b = new typename Engine::FrElement[10];
    for (int i = 0; i < 10; i++) {
        randombytes_buf((void *)&(b[i]), n8r);
        // b[i] = E.fr.random();
    }

    // TODO: remove this
    // they shouldn't really be random for testing
    for (int i = 0; i < 10; i++) {
        aux2 = E.fr.zero();
        for (int j = 0; j < 100; j++) {
            aux1 = getWitness(j+i);
            E.fr.add(aux2, aux2, aux1);
        }
        b[i] = aux2;
        LOG_DEBUG(E.fr.toString(b[i]).c_str());
    }

    typename Engine::FrElement *pol_a, *A4;
    to4T(A, domainSize, {b[2], b[1]}, pol_a, A4);
    typename Engine::FrElement *pol_b, *B4;
    to4T(B, domainSize, {b[4], b[3]}, pol_b, B4);
    typename Engine::FrElement *pol_c, *C4;
    to4T(C, domainSize, {b[6], b[5]}, pol_c, C4);

    typename Engine::G1PointAffine proof_A = expTau(pol_a, domainSize+2);
    typename Engine::G1PointAffine proof_B = expTau(pol_b, domainSize+2);
    typename Engine::G1PointAffine proof_C = expTau(pol_c, domainSize+2);

    // Round 2

    uint8_t *transcript1 = new uint8_t[64*3];
    G1toRprUncompressed(transcript1, 0, proof_A);
    G1toRprUncompressed(transcript1, 64, proof_B);
    G1toRprUncompressed(transcript1, 128, proof_C);
    // Add stuff to transcript
    /*
    LOG_DEBUG("transcript");
    LOG_DEBUG(std::to_string((int)transcript1[0]));
    LOG_DEBUG(std::to_string((int)transcript1[1]));
    LOG_DEBUG(std::to_string((int)transcript1[2]));
    LOG_DEBUG(std::to_string((int)transcript1[3]));
    LOG_DEBUG(std::to_string((int)transcript1[4]));
    */

    typename Engine::FrElement beta = hashToFr(transcript1, 64*3);

/*
    LOG_TRACE("Start Initializing a b c A");
    auto a = new typename Engine::FrElement[domainSize];
    auto b = new typename Engine::FrElement[domainSize];
    auto c = new typename Engine::FrElement[domainSize];

    #pragma omp parallel for
    for (u_int32_t i=0; i<domainSize; i++) {
        E.fr.copy(a[i], E.fr.zero());
        E.fr.copy(b[i], E.fr.zero());
    }

    LOG_TRACE("Processing coefs");
    #define NLOCKS 1024
    omp_lock_t locks[NLOCKS];
    for (int i=0; i<NLOCKS; i++) omp_init_lock(&locks[i]);
    #pragma omp parallel for 
    for (u_int64_t i=0; i<nCoefs; i++) {
        typename Engine::FrElement *ab = (coefs[i].m == 0) ? a : b;
        typename Engine::FrElement aux;

        E.fr.mul(
            aux,
            wtns[coefs[i].s],
            coefs[i].coef
        );

        omp_set_lock(&locks[coefs[i].c % NLOCKS]);
        E.fr.add(
            ab[coefs[i].c],
            ab[coefs[i].c],
            aux
        );
        omp_unset_lock(&locks[coefs[i].c % NLOCKS]);
    }
    for (int i=0; i<NLOCKS; i++) omp_destroy_lock(&locks[i]);


    LOG_TRACE("Calculating c");
    #pragma omp parallel for
    for (u_int32_t i=0; i<domainSize; i++) {
        E.fr.mul(
            c[i],
            a[i],
            b[i]
        );
    }

    LOG_TRACE("Initializing fft");
    u_int32_t domainPower = fft->log2(domainSize);

    LOG_TRACE("Start iFFT A");
    fft->ifft(a, domainSize);
    LOG_TRACE("a After ifft:");
    LOG_DEBUG(E.fr.toString(a[0]).c_str());
    LOG_DEBUG(E.fr.toString(a[1]).c_str());
    LOG_TRACE("Start Shift A");
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(a[i], a[i], fft->root(domainPower+1, i));
    }
    LOG_TRACE("a After shift:");
    LOG_DEBUG(E.fr.toString(a[0]).c_str());
    LOG_DEBUG(E.fr.toString(a[1]).c_str());
    LOG_TRACE("Start FFT A");
    fft->fft(a, domainSize);
    LOG_TRACE("a After fft:");
    LOG_DEBUG(E.fr.toString(a[0]).c_str());
    LOG_DEBUG(E.fr.toString(a[1]).c_str());
    LOG_TRACE("Start iFFT B");
    fft->ifft(b, domainSize);
    LOG_TRACE("b After ifft:");
    LOG_DEBUG(E.fr.toString(b[0]).c_str());
    LOG_DEBUG(E.fr.toString(b[1]).c_str());
    LOG_TRACE("Start Shift B");
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(b[i], b[i], fft->root(domainPower+1, i));
    }
    LOG_TRACE("b After shift:");
    LOG_DEBUG(E.fr.toString(b[0]).c_str());
    LOG_DEBUG(E.fr.toString(b[1]).c_str());
    LOG_TRACE("Start FFT B");
    fft->fft(b, domainSize);
    LOG_TRACE("b After fft:");
    LOG_DEBUG(E.fr.toString(b[0]).c_str());
    LOG_DEBUG(E.fr.toString(b[1]).c_str());

    LOG_TRACE("Start iFFT C");
    fft->ifft(c, domainSize);
    LOG_TRACE("c After ifft:");
    LOG_DEBUG(E.fr.toString(c[0]).c_str());
    LOG_DEBUG(E.fr.toString(c[1]).c_str());
    LOG_TRACE("Start Shift C");
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(c[i], c[i], fft->root(domainPower+1, i));
    }
    LOG_TRACE("c After shift:");
    LOG_DEBUG(E.fr.toString(c[0]).c_str());
    LOG_DEBUG(E.fr.toString(c[1]).c_str());
    LOG_TRACE("Start FFT C");
    fft->fft(c, domainSize);
    LOG_TRACE("c After fft:");
    LOG_DEBUG(E.fr.toString(c[0]).c_str());
    LOG_DEBUG(E.fr.toString(c[1]).c_str());

    LOG_TRACE("Start ABC");
    #pragma omp parallel for
    for (u_int64_t i=0; i<domainSize; i++) {
        E.fr.mul(a[i], a[i], b[i]);
        E.fr.sub(a[i], a[i], c[i]);
        E.fr.fromMontgomery(a[i], a[i]);
    }
    LOG_TRACE("abc:");
    LOG_DEBUG(E.fr.toString(a[0]).c_str());
    LOG_DEBUG(E.fr.toString(a[1]).c_str());

    delete b;
    delete c;

    LOG_TRACE("Start Multiexp H");
    typename Engine::G1Point pih;
    E.g1.multiMulByScalar(pih, pointsH, (uint8_t *)a, sizeof(a[0]), domainSize);
    std::ostringstream ss1;
    ss1 << "pih: " << E.g1.toString(pih);
    LOG_DEBUG(ss1);

    delete a;

    LOG_TRACE("Start Multiexp A");
    uint32_t sW = sizeof(wtns[0]);
    typename Engine::G1Point pi_a;
    E.g1.multiMulByScalar(pi_a, pointsA, (uint8_t *)wtns, sW, nVars);
    std::ostringstream ss2;
    ss2 << "pi_a: " << E.g1.toString(pi_a);
    LOG_DEBUG(ss2);

    LOG_TRACE("Start Multiexp B1");
    typename Engine::G1Point pib1;
    E.g1.multiMulByScalar(pib1, pointsB1, (uint8_t *)wtns, sW, nVars);
    std::ostringstream ss3;
    ss3 << "pib1: " << E.g1.toString(pib1);
    LOG_DEBUG(ss3);

    LOG_TRACE("Start Multiexp B2");
    typename Engine::G2Point pi_b;
    E.g2.multiMulByScalar(pi_b, pointsB2, (uint8_t *)wtns, sW, nVars);
    std::ostringstream ss4;
    ss4 << "pi_b: " << E.g2.toString(pi_b);
    LOG_DEBUG(ss4);

    LOG_TRACE("Start Multiexp C");
    typename Engine::G1Point pi_c;
    E.g1.multiMulByScalar(pi_c, pointsC, (uint8_t *)((uint64_t)wtns + (nPublic +1)*sW), sW, nVars-nPublic-1);
    std::ostringstream ss5;
    ss5 << "pi_c: " << E.g1.toString(pi_c);
    LOG_DEBUG(ss5);

    typename Engine::FrElement r;
    typename Engine::FrElement s;
    typename Engine::FrElement rs;

    E.fr.copy(r, E.fr.zero());
    E.fr.copy(s, E.fr.zero());

    randombytes_buf((void *)&(r.v[0]), sizeof(r)-1);
    randombytes_buf((void *)&(s.v[0]), sizeof(s)-1);

    typename Engine::G1Point p1;
    typename Engine::G2Point p2;

    E.g1.add(pi_a, pi_a, vk_alpha1);
    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&r, sizeof(r));
    E.g1.add(pi_a, pi_a, p1);

    E.g2.add(pi_b, pi_b, vk_beta2);
    E.g2.mulByScalar(p2, vk_delta2, (uint8_t *)&s, sizeof(s));
    E.g2.add(pi_b, pi_b, p2);

    E.g1.add(pib1, pib1, vk_beta1);
    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&s, sizeof(s));
    E.g1.add(pib1, pib1, p1);

    E.g1.add(pi_c, pi_c, pih);

    E.g1.mulByScalar(p1, pi_a, (uint8_t *)&s, sizeof(s));
    E.g1.add(pi_c, pi_c, p1);

    E.g1.mulByScalar(p1, pib1, (uint8_t *)&r, sizeof(r));
    E.g1.add(pi_c, pi_c, p1);

    E.fr.mul(rs, r, s);
    E.fr.toMontgomery(rs, rs);

    E.g1.mulByScalar(p1, vk_delta1, (uint8_t *)&rs, sizeof(rs));
    E.g1.sub(pi_c, pi_c, p1);

    */
    Proof<Engine> *p = new Proof<Engine>(Engine::engine);
    /*
    E.g1.copy(p->A, pi_a);
    E.g2.copy(p->B, pi_b);
    E.g1.copy(p->C, pi_c);
    */

    return std::unique_ptr<Proof<Engine>>(p);
}

template <typename Engine>
std::string Proof<Engine>::toJsonStr() {

    std::ostringstream ss;
    /*
    ss << "{ \"pi_a\":[\"" << E.f1.toString(A.x) << "\",\"" << E.f1.toString(A.y) << "\",\"1\"], ";
    ss << " \"pi_b\": [[\"" << E.f1.toString(B.x.a) << "\",\"" << E.f1.toString(B.x.b) << "\"],[\"" << E.f1.toString(B.y.a) << "\",\"" << E.f1.toString(B.y.b) << "\"], [\"1\",\"0\"]], ";
    ss << " \"pi_c\": [\"" << E.f1.toString(C.x) << "\",\"" << E.f1.toString(C.y) << "\",\"1\"], ";
    */
    ss << "{ \"protocol\":\"plonk\" }";
        
    return ss.str();
}

template <typename Engine>
json Proof<Engine>::toJson() {

    json p;

/*    p["pi_a"] = {};
    p["pi_a"].push_back(E.f1.toString(A.x) );
    p["pi_a"].push_back(E.f1.toString(A.y) );
    p["pi_a"].push_back("1" );


    json x2;
    x2.push_back(E.f1.toString(B.x.a));
    x2.push_back(E.f1.toString(B.x.b));
    json y2;
    y2.push_back(E.f1.toString(B.y.a));
    y2.push_back(E.f1.toString(B.y.b));
    json z2;
    z2.push_back("1");
    z2.push_back("0");
    p["pi_b"] = {};
    p["pi_b"].push_back(x2);
    p["pi_b"].push_back(y2);
    p["pi_b"].push_back(z2);

    p["pi_c"] = {};
    p["pi_c"].push_back(E.f1.toString(C.x) );
    p["pi_c"].push_back(E.f1.toString(C.y) );
    p["pi_c"].push_back("1" );
*/
    p["protocol"] = "plonk";
            
    return p;
}


} // namespace