#include "openfhe.h"

using namespace lbcrypto;

int main() {
    uint32_t multDepth = 2;
    uint32_t scaleModSize = 50;
    SecurityLevel securityLevel = HEStd_128_classic;

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(multDepth);
    parameters.SetScalingModSize(scaleModSize);
    parameters.SetBatchSize(4);
    parameters.SetSecurityLevel(securityLevel);
    parameters.SetKeySwitchTechnique(BATCHED); // new key switching technique
    
    std::cout << "Using Key Switching Technique: " << parameters.GetKeySwitchTechnique() << std::endl;
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    std::cout << "CKKS scheme is using ring dimension "
              << cc->GetRingDimension() << std::endl
              << std::endl;
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);
    cc->Enable(KEYSWITCH);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalRotateKeyGen()

    std::vector<double> msg1 = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> msg2 = {5.0, 6.0, 7.0, 8.0};

    std::cout << "msg1: " << msg1 << std::endl;
    std::cout << "msg2: " << msg2 << std::endl;
    std::cout << std::endl;

    Plaintext ptx1 = cc->MakeCKKSPackedPlaintext(msg1);
    Plaintext ptx2 = cc->MakeCKKSPackedPlaintext(msg2);

    auto ctx1 = cc->Encrypt(keyPair.publicKey, ptx1);
    auto ctx2 = cc->Encrypt(keyPair.publicKey, ptx2);
    std::cout << "Element Size: " << ctx1->GetElements().size() << std::endl;


    Plaintext decrypted_ptx1, decrypted_ptx2;
    cc->Decrypt(keyPair.secretKey, ctx1, &decrypted_ptx1);
    cc->Decrypt(keyPair.secretKey, ctx2, &decrypted_ptx2);

    decrypted_ptx1->SetLength(4);
    decrypted_ptx2->SetLength(4);

    std::vector<double> decrypted_msg1 = decrypted_ptx1->GetRealPackedValue();
    std::vector<double> decrypted_msg2 = decrypted_ptx2->GetRealPackedValue();

    std::cout << "decrypted msg1: " << decrypted_msg1 << std::endl;
    std::cout << "decrypted msg2: " << decrypted_msg2 << std::endl;
    std::cout << std::endl;

    auto ctx_mult = cc->EvalMult(ctx1, ctx2);
    std::cout << "Element Size: " << ctx_mult->GetElements().size() << std::endl;
    auto ctx_relin = cc->Relinearize(ctx_mult);
    std::cout << "Element Size: " << ctx_relin->GetElements().size() << std::endl;
 
    Plaintext ptx_mult;
    cc->Decrypt(keyPair.secretKey, ctx_relin, &ptx_mult);
    ptx_mult->SetLength(4);
    std::vector<double> msg_mult = ptx_mult->GetRealPackedValue();
    std::cout << "msg1*msg2: " << msg_mult << std::endl;
    std::cout << std::endl;

    return 0;
}