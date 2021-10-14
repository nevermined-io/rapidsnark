#include <sodium.h>

#include "logger.hpp"

using namespace CPlusPlusLogging;

namespace Plonk {


u_int32_t readUInt32(void *buf, uint32_t pos) {
    return *((u_int32_t *)((u_int64_t)buf + pos));
}


template <typename Engine>
std::string Prover<Engine>::prove(typename Engine::FrElement *_wtns) {
    // std::cerr << "brokne???\n";
    wtns = _wtns;
    // std::cerr << "brokne!\n";
    // E.fr.toString(X_2.x);
    E.f2.toString(X_2.x);
    // LOG_DEBUG(E.f2.toString(X_2.x).c_str());
    wtns[0] = E.fr.zero();
    // std::cerr << "engine???\n";

    // Calculate additions
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
        B[i] = E.fr.zero();
        C[i] = E.fr.zero();
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
    typename Engine::FrElement *ch_b = new typename Engine::FrElement[10];
    for (int i = 0; i < 10; i++) {
        ch_b[i] = E.fr.zero();
        randombytes_buf((void *)&(ch_b[i]), n8r-1);
        E.fr.mul(ch_b[i], ch_b[i], ch_b[i]);
        LOG_DEBUG(E.fr.toString(ch_b[i]).c_str());
    }

    // TODO: remove this
    // they shouldn't really be random for testing
    /*
    for (int i = 0; i < 10; i++) {
        aux2 = E.fr.zero();
        for (int j = 0; j < 100; j++) {
            aux1 = getWitness(j+i);
            E.fr.add(aux2, aux2, aux1);
        }
        // E.fr.mul(aux2, aux2, aux1);
        ch_b[i] = aux2;
        LOG_DEBUG(E.fr.toString(ch_b[i]).c_str());
    }*/

    typename Engine::FrElement *pol_a, *A4;
    to4T(A, domainSize, {ch_b[2], ch_b[1]}, pol_a, A4);
    typename Engine::FrElement *pol_b, *B4;
    to4T(B, domainSize, {ch_b[4], ch_b[3]}, pol_b, B4);
    typename Engine::FrElement *pol_c, *C4;
    to4T(C, domainSize, {ch_b[6], ch_b[5]}, pol_c, C4);

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
    LOG_DEBUG("beta");
    LOG_DEBUG(E.fr.toString(beta).c_str());

    uint8_t *transcript2 = new uint8_t[32];
    FrtoRprBE(transcript2, 0, beta);
    /*
    LOG_DEBUG("transcript");
    LOG_DEBUG(std::to_string((int)transcript2[0]));
    LOG_DEBUG(std::to_string((int)transcript2[1]));
    LOG_DEBUG(std::to_string((int)transcript2[2]));
    LOG_DEBUG(std::to_string((int)transcript2[3]));
    LOG_DEBUG(std::to_string((int)transcript2[4]));
    */

    typename Engine::FrElement gamma = hashToFr(transcript2, 32);
    LOG_DEBUG("gamma");
    LOG_DEBUG(E.fr.toString(gamma).c_str());

    typename Engine::FrElement *numArr = new typename Engine::FrElement[domainSize];
    typename Engine::FrElement *denArr = new typename Engine::FrElement[domainSize];
    typename Engine::FrElement *idenArr = new typename Engine::FrElement[domainSize];

    numArr[0] = E.fr.one();
    denArr[0] = E.fr.one();

    typename Engine::FrElement w = E.fr.one();
    typename Engine::FrElement tmp1, tmp2, tmp3, k1_mont, k2_mont;

    LOG_DEBUG("k1");
    E.fr.fromMontgomery(k1_mont, k1);
    LOG_DEBUG(E.fr.toString(k1_mont).c_str());
    LOG_DEBUG("k2");
    E.fr.fromMontgomery(k2_mont, k2);
    LOG_DEBUG(E.fr.toString(k2_mont).c_str());

    int power = log2(domainSize);

    typename Engine::FrElement *sigmaBuff = new typename Engine::FrElement[domainSize*12];

    int o = domainSize;
    for (int i = 0; i < domainSize*4; i++) {
        sigmaBuff[i] = sigmaData[o+i];
    }
    o += domainSize*5;
    for (int i = 0; i < domainSize*4; i++) {
        sigmaBuff[i+domainSize*4] = sigmaData[o+i];
    }
    o += domainSize*5;
    for (int i = 0; i < domainSize*4; i++) {
        sigmaBuff[i+domainSize*8] = sigmaData[o+i];
    }

    LOG_DEBUG("sigma buffer");
    LOG_DEBUG(E.fr.toString(sigmaBuff[0]).c_str());

    for (int i=0; i<domainSize; i++) {
    // for (int i=0; i<100; i++) {
        auto n1 = A[i];
        E.fr.mul(tmp1, beta, w);
        E.fr.add( n1, n1, tmp1 );
        E.fr.add( n1, n1, gamma );

        auto n2 = B[i];
        E.fr.mul(tmp1, beta, w);
        E.fr.mul(tmp1, tmp1, k1_mont);
        E.fr.add(n2, n2, tmp1);
        E.fr.add(n2, n2, gamma );

        auto n3 = C[i];
        E.fr.mul(tmp1, beta, w);
        E.fr.mul(tmp1, tmp1, k2_mont);
        E.fr.add(n3, n3, tmp1);
        E.fr.add(n3, n3, gamma );

        typename Engine::FrElement num;
        E.fr.mul(num, n2, n3);
        E.fr.mul(num, n1, num);

        auto d1 = A[i];
        E.fr.mul(tmp1, sigmaBuff[i*4], beta); 
        E.fr.add(d1, d1, tmp1);
        E.fr.add(d1, d1, gamma);
        // LOG_DEBUG("d1");
        // LOG_DEBUG(E.fr.toString(d1).c_str());

        auto d2 = B[i];
        E.fr.mul(tmp1, sigmaBuff[(domainSize+i)*4], beta); 
        E.fr.add(d2, d2, tmp1);
        E.fr.add(d2, d2, gamma);
        // LOG_DEBUG("d2");
        // LOG_DEBUG(E.fr.toString(d2).c_str());

        auto d3 = C[i];
        E.fr.mul(tmp1, sigmaBuff[(domainSize*2+i)*4], beta); 
        E.fr.add(d3, d3, tmp1);
        E.fr.add(d3, d3, gamma);
        // LOG_DEBUG("d3");
        // LOG_DEBUG(E.fr.toString(d3).c_str());

        typename Engine::FrElement den;
        E.fr.mul(den, d2, d3);
        E.fr.mul(den, d1, den);
        // LOG_DEBUG("den 1");
        // LOG_DEBUG(E.fr.toString(den).c_str());

        E.fr.mul(numArr[(i+1) % domainSize], numArr[i], num);
        E.fr.mul(denArr[(i+1) % domainSize], denArr[i], den);

        E.fr.mul(w, w, frw[power]);
    }

    LOG_DEBUG("final w");
    LOG_DEBUG(E.fr.toString(w).c_str());

    LOG_DEBUG("num 0");
    LOG_DEBUG(E.fr.toString(numArr[0]).c_str());
    LOG_DEBUG("den 0");
    LOG_DEBUG(E.fr.toString(denArr[0]).c_str());

    // Invert denominators
    for (int i = 0; i < domainSize; i++) {
        E.fr.inv(denArr[i], denArr[i]);
    }
    LOG_DEBUG("inv den 0");
    LOG_DEBUG(E.fr.toString(idenArr[0]).c_str());

    // E.fr.mul(tmp1, denArr[0], idenArr[0]);
    // LOG_DEBUG("inv check");
    // LOG_DEBUG(E.fr.toString(tmp1).c_str());

    for (int i=0; i<domainSize; i++) {
        E.fr.mul(numArr[i], numArr[i], denArr[i]);
    }

    auto Z = numArr;

    typename Engine::FrElement *pol_z, *Z4;
    to4T(Z, domainSize, {ch_b[9], ch_b[8], ch_b[7]}, pol_z, Z4);

    /*
    for (int i = 0; i < 10; i++) {
        std::cerr << "Z4[" << i << "] = " << E.fr.toString(Z4[i]) << "\n";
    }*/

    typename Engine::G1PointAffine proof_Z = expTau(pol_z, domainSize+3);

    // Round 3

    typename Engine::FrElement *QM4 = new typename Engine::FrElement[domainSize*4];
    for (int i = 0; i < domainSize*4; i++) {
        QM4[i] = qmData[i+domainSize];
    }
    typename Engine::FrElement *QL4 = new typename Engine::FrElement[domainSize*4];
    for (int i = 0; i < domainSize*4; i++) {
        QL4[i] = qlData[i+domainSize];
    }
    typename Engine::FrElement *QR4 = new typename Engine::FrElement[domainSize*4];
    for (int i = 0; i < domainSize*4; i++) {
        QR4[i] = qrData[i+domainSize];
    }
    typename Engine::FrElement *QO4 = new typename Engine::FrElement[domainSize*4];
    for (int i = 0; i < domainSize*4; i++) {
        QO4[i] = qoData[i+domainSize];
    }
    typename Engine::FrElement *QC4 = new typename Engine::FrElement[domainSize*4];
    for (int i = 0; i < domainSize*4; i++) {
        QC4[i] = qcData[i+domainSize];
    }
    LOG_DEBUG("QM4[0]");
    LOG_DEBUG(E.fr.toString(QM4[0]).c_str());
    LOG_DEBUG(E.fr.toString(QM4[1]).c_str());

    uint8_t *transcript3 = new uint8_t[64];
    G1toRprUncompressed(transcript3, 0, proof_Z);

    typename Engine::FrElement alpha = hashToFr(transcript3, 64);

    typename Engine::FrElement *T = new typename Engine::FrElement[domainSize*4];
    typename Engine::FrElement *Tz = new typename Engine::FrElement[domainSize*4];

    w = E.fr.one();

    for (int i=0; i<domainSize*4; i++) {
        // if ((i%4096 == 0)&&(logger)) logger.debug(`calculating t ${i}/${zkey.domainSize*4}`);

        typename Engine::FrElement a = A4[i];
        typename Engine::FrElement b = B4[i];
        typename Engine::FrElement c = C4[i];
        typename Engine::FrElement z = Z4[i];
        typename Engine::FrElement zw = Z4[(i+domainSize*4+4)%(domainSize*4)];
        /*
        if (i == 0) {
            std::cerr << "index " << (i+domainSize*4+4)%(domainSize*4) << "\n";
            std::cerr << "zw " << E.fr.toString(zw) << "\n";
        }*/
        typename Engine::FrElement qm = QM4[i];
        typename Engine::FrElement ql = QL4[i];
        typename Engine::FrElement qr = QR4[i];
        typename Engine::FrElement qo = QO4[i];
        typename Engine::FrElement qc = QC4[i];
        typename Engine::FrElement s1 = sigmaBuff[i];
        typename Engine::FrElement s2 = sigmaBuff[i+domainSize*4];
        typename Engine::FrElement s3 = sigmaBuff[i+domainSize*8];

        typename Engine::FrElement ap, bp, cp, w2, zp, wW, wW2, zWp;
        E.fr.mul(tmp1, ch_b[1], w);
        E.fr.add(ap, ch_b[2], tmp1);
        E.fr.mul(tmp1, ch_b[3], w);
        E.fr.add(bp, ch_b[4], tmp1);
        E.fr.mul(tmp1, ch_b[5], w);
        E.fr.add(cp, ch_b[6], tmp1);
        E.fr.square(w2, w);
        E.fr.mul(tmp1, ch_b[7], w2);
        E.fr.mul(tmp2, ch_b[8], w);
        E.fr.add(tmp1, tmp1, tmp2);
        E.fr.add(zp, tmp1, ch_b[9]);
        E.fr.mul(wW, w, frw[power]);
        E.fr.square(wW2, wW);
        E.fr.mul(tmp1, ch_b[7], wW2);
        E.fr.mul(tmp2, ch_b[8], wW);
        E.fr.add(tmp3, tmp1, tmp2);
        E.fr.add(zWp, tmp3, ch_b[9]);

        auto pl = E.fr.zero();
        for (int j=0; j<nPublic; j++) {
            E.fr.mul(tmp1, polData[j*5*domainSize+ domainSize+ i], A[j]);
            E.fr.sub(pl, pl, tmp1);
        }
        /*
        if (i == 0) {
            std::cerr << "pl " << E.fr.toString(pl) << "\n";
        }*/

        typename Engine::FrElement e1, e1z;
        mul2(a, b, ap, bp, i%4, e1, e1z);
        E.fr.mul(e1, e1, qm);
        E.fr.mul(e1z, e1z, qm);

        /*
        if (i == 0) {
            std::cerr << "e1 " << E.fr.toString(e1) << "\n";
            std::cerr << "e1z " << E.fr.toString(e1z) << "\n";
            std::cerr << "ql " << E.fr.toString(ql) << "\n";
            std::cerr << "a " << E.fr.toString(a) << "\n";
            std::cerr << "ap " << E.fr.toString(ap) << "\n";
        }*/
        E.fr.mul(tmp1, a, ql);
        E.fr.add(e1, e1, tmp1);
        E.fr.mul(tmp1, ap, ql);
        E.fr.add(e1z, e1z, tmp1);

        E.fr.mul(tmp1, b, qr);
        E.fr.add(e1, e1, tmp1);
        E.fr.mul(tmp1, bp, qr);
        E.fr.add(e1z, e1z, tmp1);

        E.fr.mul(tmp1, c, qo);
        E.fr.add(e1, e1, tmp1);
        E.fr.mul(tmp1, cp, qo);
        E.fr.add(e1z, e1z, tmp1);

        E.fr.add(e1, e1, pl);
        E.fr.add(e1, e1, qc);

        /*
        if (i == 0) {
            std::cerr << "e1 " << E.fr.toString(e1) << "\n";
            std::cerr << "e1z " << E.fr.toString(e1z) << "\n";
        }*/
        typename Engine::FrElement e2, e2z, betaw, e2a, e2b, e2c, e2d;
        E.fr.mul(betaw, beta, w);
        e2a =a;
        E.fr.add(e2a, e2a, betaw);
        E.fr.add(e2a, e2a, gamma);

        e2b =b;
        E.fr.mul(tmp1, betaw, k1_mont);
        E.fr.add(e2b, e2b, tmp1);
        E.fr.add(e2b, e2b, gamma);

        e2c =c;
        E.fr.mul(tmp1, betaw, k2_mont);
        E.fr.add(e2c, e2c, tmp1);
        E.fr.add(e2c, e2c, gamma);

        e2d = z;

        mul4(e2a, e2b, e2c, e2d, ap, bp, cp, zp, i%4, e2, e2z);
        /*
        if (i == 0) {
            std::cerr << "e2 " << E.fr.toString(e2) << "\n";
            std::cerr << "e2z " << E.fr.toString(e2z) << "\n";
        }*/
        E.fr.mul(e2, e2, alpha);
        E.fr.mul(e2z, e2z, alpha);

        typename Engine::FrElement e3, e3z, e3a, e3b, e3c, e3d;
        e3a = a;
        E.fr.mul(tmp1, beta, s1);
        E.fr.add(e3a, e3a, tmp1);
        E.fr.add(e3a, e3a, gamma);

        e3b = b;
        E.fr.mul(tmp1, beta, s2);
        E.fr.add(e3b, e3b,tmp1);
        E.fr.add(e3b, e3b, gamma);

        e3c = c;
        E.fr.mul(tmp1, beta, s3);
        E.fr.add(e3c, e3c, tmp1);
        E.fr.add(e3c, e3c, gamma);

        e3d = zw;
        /*
        if (i == 0) {
            std::cerr << "e3a " << E.fr.toString(e3a) << "\n";
            std::cerr << "e3b " << E.fr.toString(e3b) << "\n";
            std::cerr << "e3c " << E.fr.toString(e3c) << "\n";
            std::cerr << "e3d " << E.fr.toString(e3d) << "\n";
            std::cerr << "zWp " << E.fr.toString(zWp) << "\n";
        }*/
        mul4(e3a, e3b, e3c, e3d, ap, bp, cp, zWp, i%4, e3, e3z);
        /* if (i == 0) {
            std::cerr << "e3 " << E.fr.toString(e3) << "\n";
            std::cerr << "e3z " << E.fr.toString(e3z) << "\n";
        }*/

        E.fr.mul(e3, e3, alpha);
        E.fr.mul(e3z, e3z, alpha);

        typename Engine::FrElement e4, e4z;
        tmp1 = E.fr.one();
        E.fr.sub(e4, z, tmp1);
        E.fr.mul(e4, e4, polData[domainSize + i]);
        E.fr.mul(tmp1, alpha, alpha);
        E.fr.mul(e4, e4, tmp1);

        E.fr.mul(e4z, zp, polData[domainSize + i]);
        E.fr.mul(tmp1, alpha, alpha);
        E.fr.mul(e4z, e4z, tmp1);

        typename Engine::FrElement e, ez;
        E.fr.add(tmp1, e1, e2);
        E.fr.sub(tmp2, tmp1, e3);
        E.fr.add(e, tmp2, e4);
        E.fr.add(tmp1, e1z, e2z);
        E.fr.sub(tmp2, tmp1, e3z);
        E.fr.add(ez, tmp2, e4z);

        T[i] = e;
        Tz[i] = ez;

        E.fr.mul(w, w, frw[power+2]);
    }

    /*
    std::cerr << "final w " << E.fr.toString(w) << "\n";
    std::cerr << "T[1001] " << E.fr.toString(T[1001]) << "\n";
    std::cerr << "Tz[1001] " << E.fr.toString(Tz[1001]) << "\n";
    */

    fft->ifft(T, domainSize*4);
    typename Engine::FrElement *t = T;
    // std::cerr << "t[1001] " << E.fr.toString(t[1001]) << "\n";

    // dividing T/Z    
    for (int i=0; i<domainSize; i++) {
        E.fr.neg(t[i], t[i]);
    }

    for (int i=domainSize; i<domainSize*4; i++) {
        E.fr.sub(t[i], t[i-domainSize], t[i]);
        /*
        if (i > (zkey.domainSize*3 -4) ) {
            if (!Fr.isZero(a)) {
                throw new Error("T Polynomial is not divisible");
            }
        }*/
    }
    // std::cerr << "t[1001+size] " << E.fr.toString(t[1001+domainSize]) << "\n";

    fft->ifft(Tz, domainSize*4);
    typename Engine::FrElement *tz = Tz;
    for (int i=0; i<domainSize*4; i++) {
        // auto a = tz[i];
        if (i > domainSize*3 +5) {
            /*
            if (!Fr.isZero(a)) {
                throw new Error("Tz Polynomial is not well calculated");
            }*/
        } else {
            E.fr.add(t[i], t[i], tz[i]);
        }
    }
    // std::cerr << "t[1001+size] " << E.fr.toString(t[1001+domainSize]) << "\n";

    typename Engine::FrElement *pol_t = t; // size: domainSize*3+6

    auto proof_T1 = expTau(t, domainSize);
    auto proof_T2 = expTau(t+domainSize, domainSize);
    auto proof_T3 = expTau(t+domainSize*2, domainSize+6);

    // Round 4
    typename Engine::FrElement *pol_qm = qmData;
    typename Engine::FrElement *pol_ql = qlData;
    typename Engine::FrElement *pol_qr = qrData;
    typename Engine::FrElement *pol_qo = qoData;
    typename Engine::FrElement *pol_qc = qcData;
    typename Engine::FrElement *pol_s1 = sigmaData + 0*domainSize;
    typename Engine::FrElement *pol_s2 = sigmaData + 5*domainSize;
    typename Engine::FrElement *pol_s3 = sigmaData + 10*domainSize;

    uint8_t *transcript4 = new uint8_t[64*3];
    G1toRprUncompressed(transcript4, 0, proof_T1);
    G1toRprUncompressed(transcript4, 64, proof_T2);
    G1toRprUncompressed(transcript4, 128, proof_T3);
    auto xi = hashToFr(transcript4, 64*3);

    auto proof_eval_a = evalPol(pol_a, xi, domainSize+2);
    auto proof_eval_b = evalPol(pol_b, xi, domainSize+2);
    auto proof_eval_c = evalPol(pol_c, xi, domainSize+2);
    auto proof_eval_s1 = evalPol(pol_s1, xi, domainSize);
    auto proof_eval_s2 = evalPol(pol_s2, xi, domainSize);
    auto proof_eval_t = evalPol(pol_t, xi, domainSize*3+6);
    E.fr.mul(tmp1, xi, frw[power]);
    auto proof_eval_zw = evalPol(pol_z, tmp1, domainSize+3);

    typename Engine::FrElement betaxi, coef_ab;
    E.fr.mul(coef_ab, proof_eval_a, proof_eval_b);
    
    E.fr.mul(betaxi, beta, xi);

    auto e2a = proof_eval_a;
    E.fr.add(e2a, e2a, betaxi);
    E.fr.add(e2a, e2a, gamma);

    auto e2b = proof_eval_b;
    E.fr.mul(tmp1, betaxi, k1_mont);
    E.fr.add(e2b, e2b, tmp1);
    E.fr.add(e2b, e2b, gamma);

    auto e2c = proof_eval_c;
    E.fr.mul(tmp1, betaxi, k2_mont);
    E.fr.add(e2c, e2c, tmp1);
    E.fr.add(e2c, e2c, gamma);

    typename Engine::FrElement e2;
    E.fr.mul(tmp1, e2a, e2b);
    E.fr.mul(tmp1, tmp1, e2c);
    E.fr.mul(e2, tmp1, alpha);

    auto e3a = proof_eval_a;
    E.fr.mul(tmp1, beta, proof_eval_s1);
    E.fr.add(e3a, e3a, tmp1);
    E.fr.add(e3a, e3a, gamma);

    auto e3b = proof_eval_b;
    E.fr.mul(tmp1, beta, proof_eval_s2);
    E.fr.add(e3b, e3b, tmp1);
    E.fr.add(e3b, e3b, gamma);

    typename Engine::FrElement e3;
    E.fr.mul(e3, e3a, e3b);
    E.fr.mul(e3, e3, beta);
    E.fr.mul(e3, e3, proof_eval_zw);
    E.fr.mul(e3, e3, alpha);

    auto xim = xi;
    for (int i=0; i<power; i++) E.fr.mul(xim, xim, xim);
    typename Engine::FrElement eval_l1, e4, coefz, coefs3, size;
    auto one = E.fr.one();
    E.fr.sub(tmp1, xi, one);
    E.fr.sub(tmp2, xim, one);
    E.fr.fromUI(size, domainSize);
    E.fr.mul(tmp3, tmp1, size);
    E.fr.div(eval_l1, tmp2, tmp3);

    E.fr.mul(tmp1, alpha, alpha);
    E.fr.mul(e4, eval_l1, tmp1);

    coefs3 = e3;
    E.fr.add(coefz, e2, e4);

    typename Engine::FrElement *pol_r = new typename Engine::FrElement[domainSize*3];

    typename Engine::FrElement v;
    for (int i = 0; i<domainSize+3; i++) {
        E.fr.mul(v, coefz, pol_z[i]);
        if (i<domainSize) {
            E.fr.mul(tmp1, coef_ab, pol_qm[i]);
            E.fr.add(v, v, tmp1);
            E.fr.mul(tmp1, proof_eval_a, pol_ql[i]);
            E.fr.add(v, v, tmp1);
            E.fr.mul(tmp1, proof_eval_b, pol_qr[i]);
            E.fr.add(v, v, tmp1);
            E.fr.mul(tmp1, proof_eval_c, pol_qo[i]);
            E.fr.add(v, v, tmp1);
            E.fr.add(v, v, pol_qc[i]);
            E.fr.mul(tmp1, coefs3, pol_s3[i]);
            E.fr.sub(v, v, tmp1);
        }
        pol_r[i] = v;
    }

    auto proof_eval_r = evalPol(pol_r, xi, domainSize+3);

    // Round 5
    int n8r = 32;
    uint8_t *transcript5 = new uint8_t[32*7];
    FrtoRprBE(transcript5, 0, proof_eval_a);
    FrtoRprBE(transcript5, n8r, proof_eval_b);
    FrtoRprBE(transcript5, n8r*2, proof_eval_c);
    FrtoRprBE(transcript5, n8r*3, proof_eval_s1);
    FrtoRprBE(transcript5, n8r*4, proof_eval_s2);
    FrtoRprBE(transcript5, n8r*5, proof_eval_zw);
    FrtoRprBE(transcript5, n8r*6, proof_eval_r);
    typename Engine::FrElement *ch_v = new typename Engine::FrElement[7];
    ch_v[1] = hashToFr(transcript5, 32*7);

    for (int i=2; i<=6; i++ ) {
        E.fr.mul(ch_v[i], ch_v[i-1], ch_v[1]);
    }

    typename Engine::FrElement *pol_wxi = new typename Engine::FrElement[domainSize+6];
    typename Engine::FrElement xi2m;

    E.fr.mul(xi2m, xim, xim);

    for (int i=0; i<domainSize+6; i++) {
        w = E.fr.zero();
        E.fr.mul(tmp1, xi2m, pol_t[domainSize*2+i]);
        E.fr.add(w, w, tmp1);

        if (i<domainSize+3) {
            E.fr.mul(tmp1, ch_v[1],  pol_r[i]);
            E.fr.add(w, w, tmp1);
        }

        if (i<domainSize+2) {
            E.fr.mul(tmp1, ch_v[2],  pol_a[i]);
            E.fr.add(w, w, tmp1);
            E.fr.mul(tmp1, ch_v[3],  pol_b[i]);
            E.fr.add(w, w, tmp1);
            E.fr.mul(tmp1, ch_v[4],  pol_c[i]);
            E.fr.add(w, w, tmp1);
        }
        
        if (i<domainSize) {
            E.fr.add(w, w, pol_t[i]);
            E.fr.mul(tmp1, xim,  pol_t[domainSize+i]);
            E.fr.add(w, w, tmp1);
            E.fr.mul(tmp1, ch_v[5],  pol_s1[i]);
            E.fr.add(w, w, tmp1);
            E.fr.mul(tmp1, ch_v[6],  pol_s2[i]);
            E.fr.add(w, w, tmp1);
        }

        pol_wxi[i] = w;
    }

    auto w0 = pol_wxi[0];
    E.fr.sub(w0, w0, proof_eval_t);
    E.fr.mul(tmp1, ch_v[1], proof_eval_r);
    E.fr.sub(w0, w0, tmp1);
    E.fr.mul(tmp1, ch_v[2], proof_eval_a);
    E.fr.sub(w0, w0, tmp1);
    E.fr.mul(tmp1, ch_v[3], proof_eval_b);
    E.fr.sub(w0, w0, tmp1);
    E.fr.mul(tmp1, ch_v[4], proof_eval_c);
    E.fr.sub(w0, w0, tmp1);
    E.fr.mul(tmp1, ch_v[5], proof_eval_s1);
    E.fr.sub(w0, w0, tmp1);
    E.fr.mul(tmp1, ch_v[6], proof_eval_s2);
    E.fr.sub(w0, w0, tmp1);
    pol_wxi[0] = w0;

    typename Engine::FrElement *pol_wxi_;
    divPol1(pol_wxi, xi, domainSize+6, pol_wxi_);

    auto proof_Wxi = expTau(pol_wxi_, domainSize+6);

    typename Engine::FrElement *pol_wxiw = new typename Engine::FrElement[domainSize+3];
    for (int i=0; i<domainSize+3; i++) {
        pol_wxiw[i] = pol_z[i];
    }
    w0 = pol_wxiw[0];
    E.fr.sub(w0, w0, proof_eval_zw);
    pol_wxiw[0] = w0;

    typename Engine::FrElement *pol_wxiw_;
    E.fr.mul(tmp1, xi, frw[power]);
    divPol1(pol_wxiw, tmp1, domainSize+3, pol_wxiw_);
    auto proof_Wxiw = expTau(pol_wxiw_, domainSize+3);

    uint8_t *res = new uint8_t[n8r*25];
    G1toRprUncompressed(res, n8r*0, proof_A);
    G1toRprUncompressed(res, n8r*2, proof_B);
    G1toRprUncompressed(res, n8r*4, proof_C);
    G1toRprUncompressed(res, n8r*6, proof_Z);
    G1toRprUncompressed(res, n8r*8, proof_T1);
    G1toRprUncompressed(res, n8r*10, proof_T2);
    G1toRprUncompressed(res, n8r*12, proof_T3);
    G1toRprUncompressed(res, n8r*14, proof_Wxi);
    G1toRprUncompressed(res, n8r*16, proof_Wxiw);
    FrtoRprBE(res, n8r*18, proof_eval_a);
    FrtoRprBE(res, n8r*19, proof_eval_b);
    FrtoRprBE(res, n8r*20, proof_eval_c);
    FrtoRprBE(res, n8r*21, proof_eval_s1);
    FrtoRprBE(res, n8r*22, proof_eval_s2);
    FrtoRprBE(res, n8r*23, proof_eval_zw);
    FrtoRprBE(res, n8r*24, proof_eval_r);
    std::stringstream ss;
    ss << "0x";
    for (int i = 0; i < n8r*25; i++) {
        if (res[i] < 16) {
            ss << "0" << std::hex << int(res[i]);
        } else {
            ss << std::hex << int(res[i]);
        }
    }
    std::string result;
    ss >> result;
    // std::cerr << result << "\n";

    return result;
}

} // namespace