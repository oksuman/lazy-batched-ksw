
#pragma once

#include <openfhe.h>
#include <memory>
#include <vector>

using namespace lbcrypto;

class MATMULT_AS24{
    public: 
        MATMULT_AS24(int dim) {
            this->d = dim;

            max_batch = this->m_cc->GetRingDimension() / 2;
            s = std::min(max_batch / d / d, d);
            B = d / s;
        } 

    private: 
        CryptoContext<DCRTPoly> m_cc;
        PublicKey<DCRTPoly> m_PublicKey;
        
        int d; //matrix dimension
        int max_batch;
        int B;
        int s;


};