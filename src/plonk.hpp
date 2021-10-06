#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "binfile_utils.hpp"
#include "fft.hpp"

extern "C" {
    #include <sha3.h>
}

#include "logger.hpp"
using namespace CPlusPlusLogging;

namespace Plonk {

    template <typename Engine>
    class Proof {
        Engine &E;
    public:
        typename Engine::G1PointAffine A;
        typename Engine::G2PointAffine B;
        typename Engine::G1PointAffine C;

        Proof(Engine &_E) : E(_E) { };
        std::string toJsonStr();
        json toJson();
    };


 #pragma pack(push, 1)
    template <typename Engine>
    struct Coef {
        u_int32_t m;
        u_int32_t c;
        u_int32_t s;
        typename Engine::FrElement coef;
    };
#pragma pack(pop)

    template <typename Engine>
    class Prover {

        Engine &E;
        u_int32_t nVars;
        u_int32_t nPublic;
        u_int32_t domainSize;
        u_int32_t nAdditions;
        u_int32_t nConstrain;
        u_int64_t n8r;
        u_int64_t n8q;
        typename Engine::G1PointAffine &Qm;
        typename Engine::G1PointAffine &Ql;
        typename Engine::G1PointAffine &Qr;
        typename Engine::G1PointAffine &Qo;
        typename Engine::G1PointAffine &Qc;
        typename Engine::G1PointAffine &S1;
        typename Engine::G1PointAffine &S2;
        typename Engine::G1PointAffine &S3;
        typename Engine::G2PointAffine &X_2;

        typename Engine::G1PointAffine *PTau;

        typename Engine::FrElement *internalWitness;
        typename Engine::FrElement *wtns;
        typename Engine::FrElement k1;
        typename Engine::FrElement k2;

        uint8_t *additions;
        uint8_t *Amap;
        uint8_t *Bmap;
        uint8_t *Cmap;

        typename Engine::FrElement *sigmaData;

        typename Engine::FrElement *qmData;
        typename Engine::FrElement *qlData;
        typename Engine::FrElement *qrData;
        typename Engine::FrElement *qoData;
        typename Engine::FrElement *qcData;
        typename Engine::FrElement *polData;
        typename Engine::FrElement frw[29];

        typename Engine::FrElement Z1[4];
        typename Engine::FrElement Z2[4];
        typename Engine::FrElement Z3[4];

        FFT<typename Engine::Fr> *fft;
    public:
        Prover(
            Engine &_E, 
            u_int32_t _nVars, 
            u_int32_t _nPublic, 
            u_int32_t _domainSize, 
            u_int32_t _nAdditions,
            u_int32_t _nConstrain,
            typename Engine::G1PointAffine &_Qm,
            typename Engine::G1PointAffine &_Ql,
            typename Engine::G1PointAffine &_Qr,
            typename Engine::G1PointAffine &_Qo,
            typename Engine::G1PointAffine &_Qc,
            typename Engine::G1PointAffine &_S1,
            typename Engine::G1PointAffine &_S2,
            typename Engine::G1PointAffine &_S3,
            typename Engine::G2PointAffine &_X_2,
            u_int32_t _n8r,
            u_int32_t _n8q,
            uint8_t *_additions,
            uint8_t *_aMap,
            uint8_t *_bMap,
            uint8_t *_cMap,
            typename Engine::G1PointAffine *_PTau,
            typename Engine::FrElement _k1,
            typename Engine::FrElement _k2,
            typename Engine::FrElement *_sigmaData,
            typename Engine::FrElement *_qmData,
            typename Engine::FrElement *_qlData,
            typename Engine::FrElement *_qrData,
            typename Engine::FrElement *_qoData,
            typename Engine::FrElement *_qcData,
            typename Engine::FrElement *_polData
        ) : 
            E(_E), 
            nVars(_nVars),
            nPublic(_nPublic),
            domainSize(_domainSize),
            nAdditions(_nAdditions),
            nConstrain(_nConstrain),
            n8r(_n8r),
            n8q(_n8q),
            Qm(_Qm),
            Ql(_Ql),
            Qr(_Qr),
            Qo(_Qo),
            Qc(_Qc),
            S1(_S1),
            S2(_S2),
            S3(_S3),
            X_2(_X_2),
            PTau(_PTau),
            additions(_additions),
            Amap(_aMap),
            Bmap(_bMap),
            Cmap(_cMap),
            k1(_k1),
            k2(_k2),
            sigmaData(_sigmaData),
            qmData(_qmData),
            qlData(_qlData),
            qrData(_qrData),
            qoData(_qoData),
            qcData(_qcData),
            polData(_polData)
        {
            fft = new FFT<typename Engine::Fr>(domainSize*4);
            internalWitness = new typename Engine::FrElement[nAdditions];
            // LOG_DEBUG("Fr.w");
            E.fr.fromString(frw[28], "19103219067921713944291392827692070036145651957329286315305642004821462161904");
            for (int i = 27; i >= 0; i--) {
                E.fr.square(frw[i], frw[i+1]);
                // LOG_DEBUG(E.fr.toString(frw[i]).c_str());
            }

            typename Engine::FrElement minus_1;
            typename Engine::FrElement minus_2;
            typename Engine::FrElement zero = E.fr.zero();
            typename Engine::FrElement one = E.fr.one();
            typename Engine::FrElement two;
            typename Engine::FrElement tmp;
            E.fr.sub(minus_1, zero, one);
            E.fr.sub(minus_2, minus_1, one);

            Z1[0] = E.fr.zero();
            E.fr.add(Z1[1], minus_1, frw[2]);
            Z1[2] = minus_2;
            E.fr.sub(Z1[3], minus_1, frw[2]);

            Z2[0] = E.fr.zero();
            E.fr.mul(tmp, minus_2, frw[2]);
            Z2[1] = tmp;
            E.fr.fromString(Z2[2], "4");
            E.fr.sub(Z2[3], zero, tmp);

            E.fr.fromString(two, "2");

            E.fr.mul(tmp, two, frw[2]);
            Z3[0] = E.fr.zero();
            E.fr.add(Z3[1], two, tmp);
            E.fr.fromString(Z3[2], "-8");
            E.fr.sub(Z3[3], two, tmp);

            LOG_DEBUG("Z3");
            LOG_DEBUG(E.fr.toString(Z3[0]).c_str());
            LOG_DEBUG(E.fr.toString(Z3[1]).c_str());
            LOG_DEBUG(E.fr.toString(Z3[2]).c_str());
            LOG_DEBUG(E.fr.toString(Z3[3]).c_str());
        };

        ~Prover() {
            delete fft;
        }

        std::unique_ptr<Proof<Engine>> prove(typename Engine::FrElement *wtns);

        void mul2(typename Engine::FrElement &a,
                typename Engine::FrElement &b,
                typename Engine::FrElement &ap,
                typename Engine::FrElement &bp,
                int p,
                typename Engine::FrElement &r,
                typename Engine::FrElement &rz
                ) {
            
            typename Engine::FrElement a_b, a_bp, ap_b, ap_bp, a0, a1;
            E.fr.mul(a_b,a,b);
            E.fr.mul(a_bp,a,bp);
            E.fr.mul(ap_b,ap,b);
            E.fr.mul(ap_bp,ap,bp);

            r = a_b;

            E.fr.add(a0, a_bp, ap_b);

            a1 = ap_bp;

            rz = a0;
            if (p > 0) {
                E.fr.mul(a0, Z1[p], a1);
                E.fr.add(rz, rz, a0);
            }
        }

        void mul4(typename Engine::FrElement &a,
                typename Engine::FrElement &b,
                typename Engine::FrElement &c,
                typename Engine::FrElement &d,
                typename Engine::FrElement &ap,
                typename Engine::FrElement &bp,
                typename Engine::FrElement &cp,
                typename Engine::FrElement &dp,
                int p,
                typename Engine::FrElement &r,
                typename Engine::FrElement &rz
                ) {

            
            typename Engine::FrElement a_b, a_bp, ap_b, ap_bp, a0, a1, a2, a3;
            typename Engine::FrElement c_d, c_dp, cp_d, cp_dp, tmp;
            E.fr.mul(a_b,a,b);
            E.fr.mul(a_bp,a,bp);
            E.fr.mul(ap_b,ap,b);
            E.fr.mul(ap_bp,ap,bp);

            E.fr.mul(c_d,c,d);
            E.fr.mul(c_dp,c,dp);
            E.fr.mul(cp_d,cp,d);
            E.fr.mul(cp_dp,cp,dp);

            E.fr.mul(r, a_b, c_d);

            E.fr.mul(a0, ap_b, c_d);
            E.fr.mul(tmp, a_bp, c_d);
            E.fr.add(a0, a0, tmp);
            E.fr.mul(tmp, a_b, cp_d);
            E.fr.add(a0, a0, tmp);
            E.fr.mul(tmp, a_b, c_dp);
            E.fr.add(a0, a0, tmp);

            E.fr.mul(a1, ap_bp, c_d);
            E.fr.mul(tmp,ap_b, cp_d);
            E.fr.add(a1, a1, tmp);
            E.fr.mul(tmp,ap_b, c_dp);
            E.fr.add(a1, a1, tmp);
            E.fr.mul(tmp,a_bp, cp_d);
            E.fr.add(a1, a1, tmp);
            E.fr.mul(tmp,a_bp, c_dp);
            E.fr.add(a1, a1, tmp);
            E.fr.mul(tmp,a_b, cp_dp);
            E.fr.add(a1, a1, tmp);

            E.fr.mul(a2, a_bp, cp_dp);
            E.fr.mul(tmp,ap_b, cp_dp);
            E.fr.add(a2, a2, tmp);
            E.fr.mul(tmp,ap_bp, c_dp);
            E.fr.add(a2, a2, tmp);
            E.fr.mul(tmp,ap_bp, cp_d);
            E.fr.add(a2, a2, tmp);

            E.fr.mul(a3, ap_bp, cp_dp);

            rz = a0;
            if (p > 0) {
                E.fr.mul(tmp, Z1[p], a1);
                E.fr.add(rz, rz, tmp);
                E.fr.mul(tmp, Z2[p], a2);
                E.fr.add(rz, rz, tmp);
                E.fr.mul(tmp, Z3[p], a3);
                E.fr.add(rz, rz, tmp);
            }

        }

        typename Engine::FrElement hashToFr(uint8_t *transcript, int size) {
            uint8_t out[32];
            sha3_HashBuffer(256, SHA3_FLAGS_KECCAK, transcript, size, out, sizeof(out));
            /*
            LOG_DEBUG("hash keccak");
            LOG_DEBUG(std::to_string((int)out[0]));
            LOG_DEBUG(std::to_string((int)out[1]));
            LOG_DEBUG(std::to_string((int)out[2]));
            LOG_DEBUG(std::to_string((int)out[3]));
            LOG_DEBUG(std::to_string((int)out[4]));
            */
            std::string tmp = "";
            typename Engine::FrElement res = E.fr.zero();
            typename Engine::FrElement tmp_num;
            typename Engine::FrElement e256;
            E.fr.fromString(e256, "256");
            for (int i = 0; i < 32; i++) {
                int j = i;
                char d1 = out[j] % 10;
                char d2 = out[j] / 10 % 10;
                char d3 = out[j] / 100 % 10;
                // std::cerr << (int)d3 << " " << (int)d2 << " " << (int)d1 << " -- " << (int)out[j] << "\n";
                tmp = "";
                tmp += ('0' + d3);
                tmp += ('0'+d2);
                tmp += ('0'+d1);
                E.fr.fromString(tmp_num, tmp);
                // std::cerr << "tmp " << E.fr.toString(tmp_num) << " from " << tmp << "\n";
                E.fr.mul(res, res, e256);
                // std::cerr << "acc " << E.fr.toString(res) << " mul " << E.fr.toString(e256) << "\n";
                E.fr.add(res, res, tmp_num);
                // std::cerr << "acc " << E.fr.toString(res) << "\n";
            }
            LOG_DEBUG("hash keccak");
            LOG_DEBUG(E.fr.toString(res).c_str());
            return res;
        }

        void G1toRprUncompressed(uint8_t *transcript1, int offset, typename Engine::G1PointAffine &p) {
            uint8_t *tmp;
            typename Engine::F1Element bm;
            E.f1.fromMontgomery(bm, p.x);
            tmp = (uint8_t*)(&bm);
            for (int i = 0; i < 32; i++) {
                transcript1[i+offset] = tmp[31-i];
            }
            E.f1.fromMontgomery(bm, p.y);
            tmp = (uint8_t*)(&bm);
            for (int i = 0; i < 32; i++) {
                transcript1[i+offset+32] = tmp[31-i];
            }

        }

        void FrtoRprBE(uint8_t *transcript1, int offset, typename Engine::FrElement &p) {
            typename Engine::FrElement bm;
            E.fr.fromMontgomery(bm, p);
            uint8_t *tmp = (uint8_t*)(&bm);
            for (int i = 0; i < 32; i++) {
                transcript1[i+offset] = tmp[31-i];
            }
        }

        void divPol1(typename Engine::FrElement *P, typename Engine::FrElement &d, int n, typename Engine::FrElement *&res) {
            res = new typename Engine::FrElement[n];
            res[n-1] = E.fr.zero();
            res[n-2] = P[n-1];
            for (int i=n-3; i>=0; i--) {
                E.fr.mul(res[i], d, res[i+1]);
                E.fr.add(res[i], res[i], P[i+1]);
            }
            /*
            if (!Fr.eq(
                P.slice(0, n8r),
                Fr.mul(
                    Fr.neg(d),
                    res.slice(0, n8r)
                )
            )) {
                throw new Error("Polinomial does not divide");
            }
            */
        }

        typename Engine::FrElement evalPol(typename Engine::FrElement *P, typename Engine::FrElement &x, int n) {
            if (n == 0) return E.fr.zero();
            auto res = P[n-1];
            for (int i=n-2; i>=0; i--) {
                E.fr.mul(res, res, x);
                E.fr.add(res, res, P[i]);
            }
            std::cerr << "evalpol " << n << " res " << E.fr.toString(res) << "\n";
            return res;
        }


        typename Engine::G1PointAffine expTau(typename Engine::FrElement *b, int size) {
            typename Engine::FrElement *bm = new typename Engine::FrElement[size];
            for (int i = 0; i < size; i++) {
                E.fr.fromMontgomery(bm[i], b[i]);
            }
            LOG_DEBUG("from montgomery");
            LOG_DEBUG(E.fr.toString(bm[0]).c_str());
            typename Engine::G1Point res;
            typename Engine::G1PointAffine res2;
            E.g1.multiMulByScalar(res, PTau, (uint8_t *)bm, sizeof(bm[0]), size);
            E.g1.copy(res2, res);
            /*
            LOG_DEBUG("size");
            LOG_DEBUG(std::to_string(size));
            */
            LOG_DEBUG("exptau result");
            LOG_DEBUG(E.f1.toString(res2.x).c_str());
            /*
            LOG_DEBUG("ptau[0]");
            LOG_DEBUG(E.f1.toString(PTau[0].x).c_str());
            LOG_DEBUG(E.f1.toString(PTau[0].y).c_str());
            LOG_DEBUG(E.f1.toString(PTau[1].x).c_str());
            LOG_DEBUG(E.f1.toString(PTau[1].y).c_str());
            */
            return res2;
        }

        void to4T(typename Engine::FrElement *A, uint32_t size, std::vector<typename Engine::FrElement> pz, typename Engine::FrElement* & a1, typename Engine::FrElement* &A4) {
            typename Engine::FrElement *a = new typename Engine::FrElement[2*size];
            
            for (int i = 0; i < size*2; i++) {
                a[i] = E.fr.zero();
            }
            for (int i = 0; i < size; i++) {
                a[i] = A[i];
                // LOG_DEBUG("hmm");
                // LOG_DEBUG(E.fr.toString(a[i]).c_str());
            }
            /*
            LOG_DEBUG("before ifft");
            LOG_DEBUG(E.fr.toString(a[0]).c_str());
            LOG_DEBUG("size");
            LOG_DEBUG(std::to_string(size));
            */
            u_int32_t domainPower = fft->log2(size);
            fft->ifft(a, size);
            LOG_DEBUG("after ifft");
            LOG_DEBUG(E.fr.toString(a[0]).c_str());

            typename Engine::FrElement *a4 = new typename Engine::FrElement[4*size];
            for (int i = 0; i < size*4; i++) {
                a4[i] = E.fr.zero();
            }
            for (int i = 0; i < size*2; i++) {
                a4[i] = a[i];
            }

            a1 = new typename Engine::FrElement[size+pz.size()];
            for (int i = 0; i < size+pz.size(); i++) {
                a1[i] = a[i];
            }

            typename Engine::FrElement tmp;
            for (int i = 0; i < pz.size(); i++) {
                E.fr.add(a1[size+i], a1[size+i], pz[i]);
                E.fr.sub(a1[i], a1[i], pz[i]);
                LOG_DEBUG("randoms");
                LOG_DEBUG(E.fr.toString(a1[i]).c_str());
                LOG_DEBUG(E.fr.toString(a1[size+i]).c_str());
            }

            /*
            LOG_DEBUG("nVars");
            LOG_DEBUG(std::to_string(nVars));
            LOG_DEBUG("nAdditions");
            LOG_DEBUG(std::to_string(nAdditions));
            */
            fft->fft(a4, size*4);
            A4 = a4;
            LOG_DEBUG("A4");
            LOG_DEBUG(E.fr.toString(A4[1]).c_str());
        }

        typename Engine::FrElement& getWitness(uint32_t idx) {
            if (idx < nVars-nAdditions) {
                return *(typename Engine::FrElement*)(wtns+idx);
            } else if (idx < nVars) {
                return *(typename Engine::FrElement*)(internalWitness + (idx - (nVars-nAdditions)));
            } else {
                return E.fr.zero();
            }
        }
    };

    template <typename Engine>
    std::unique_ptr<Prover<Engine>> makeProver(
        u_int32_t nVars, 
        u_int32_t nPublic, 
        u_int32_t domainSize,
        u_int32_t nAdditions,
        u_int32_t nConstrain,
        void *Qm,
        void *Ql,
        void *Qr,
        void *Qo,
        void *Qc,
        void *S1,
        void *S2,
        void *S3,
        void *X_2,
        u_int32_t n8r,
        u_int32_t n8q,
        void *additions,
        void *amap,
        void *bmap,
        void *cmap,
        void *PTau,
        mpz_t _k1,
        mpz_t _k2,
        void *sigmaData,
        void *qmData,
        void *qlData,
        void *qrData,
        void *qoData,
        void *qcData,
        void *polData
    ) {
        typename Engine::FrElement k1;
        typename Engine::FrElement k2;
        Engine::engine.fr.fromMpz(k1, _k1);
        Engine::engine.fr.fromMpz(k2, _k2);
        Prover<Engine> *p = new Prover<Engine>(
            Engine::engine, 
            nVars, 
            nPublic, 
            domainSize, 
            nAdditions,
            nConstrain,
            *(typename Engine::G1PointAffine *)Qm,
            *(typename Engine::G1PointAffine *)Ql,
            *(typename Engine::G1PointAffine *)Qr,
            *(typename Engine::G1PointAffine *)Qo,
            *(typename Engine::G1PointAffine *)Qc,
            *(typename Engine::G1PointAffine *)S1,
            *(typename Engine::G1PointAffine *)S2,
            *(typename Engine::G1PointAffine *)S3,
            *(typename Engine::G2PointAffine *)X_2,
            n8r,
            n8q,
            (uint8_t*)additions,
            (uint8_t*)amap,
            (uint8_t*)bmap,
            (uint8_t*)cmap,
            (typename Engine::G1PointAffine *)PTau,
            k1,
            k2,
            (typename Engine::FrElement *)sigmaData,
            (typename Engine::FrElement *)qmData,
            (typename Engine::FrElement *)qlData,
            (typename Engine::FrElement *)qrData,
            (typename Engine::FrElement *)qoData,
            (typename Engine::FrElement *)qcData,
            (typename Engine::FrElement *)polData
        );
        return std::unique_ptr< Prover<Engine> >(p);
    }
};


#include "plonk.cpp"
