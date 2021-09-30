#include <string>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "binfile_utils.hpp"
#include "fft.hpp"
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

        typename Engine::FrElement *internalWitness;
        typename Engine::FrElement *wtns;

        uint8_t *additions;
        uint8_t *Amap;
        uint8_t *Bmap;
        uint8_t *Cmap;

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
            uint8_t *_cMap
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
            additions(_additions),
            Amap(_aMap),
            Bmap(_bMap),
            Cmap(_cMap)
        {
            fft = new FFT<typename Engine::Fr>(domainSize*4);
            internalWitness = new typename Engine::FrElement[nAdditions];
        };

        ~Prover() {
            delete fft;
        }

        std::unique_ptr<Proof<Engine>> prove(typename Engine::FrElement *wtns);

        void to4T(typename Engine::FrElement *A, uint32_t size, std::vector<typename Engine::FrElement> pz, typename Engine::FrElement* & a1, typename Engine::FrElement* &A4) {
            typename Engine::FrElement *a = new typename Engine::FrElement[2*size];
            for (int i = 0; i < size*2; i++) {
                a[i] = E.fr.zero();
            }
            a[0] = A[0];
            a[1] = A[1];
            a[2] = A[1];
            // size = 16;
            for (int i = 0; i < size; i++) {
                a[i] = A[i];
                LOG_DEBUG("hmm");
                LOG_DEBUG(E.fr.toString(a[i]).c_str());
            }
            LOG_DEBUG("before ifft");
            LOG_DEBUG(E.fr.toString(a[0]).c_str());
            LOG_DEBUG("size");
            LOG_DEBUG(std::to_string(size));
            u_int32_t domainPower = fft->log2(size);
            fft->ifft(a, size);
            LOG_DEBUG("after ifft");
            LOG_DEBUG(E.fr.toString(a[0]).c_str());
            LOG_DEBUG("nVars");
            LOG_DEBUG(std::to_string(nVars));
            LOG_DEBUG("nAdditions");
            LOG_DEBUG(std::to_string(nAdditions));
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
        void *cmap
    ) {
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
            (uint8_t*)cmap
        );
        return std::unique_ptr< Prover<Engine> >(p);
    }
};


#include "plonk.cpp"
