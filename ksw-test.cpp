#include "openfhe.h"

using namespace lbcrypto;

void CKKSTest(){
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
    cc->EvalLazyRotateKeyGen(keyPair.secretKey, {1,2,3,6});

    std::vector<double> msg1 = {1.0, 2.0, 3.0, 4.0};
    std::vector<double> msg2 = {5.0, 6.0, 7.0, 8.0};
    std::vector<double> msg3 = {9.0, 10.0, 11.0, 12.0};

    std::vector<double> msk1 = {1.0, 0.0, 0.0, 0.0};
    std::vector<double> msk2 = {0.0, 1.0, 0.0, 0.0};
    std::vector<double> msk3 = {0.0, 0.0, 1.0, 0.0};

    std::cout << "msg1: " << msg1 << std::endl;
    std::cout << "msg2: " << msg2 << std::endl;
    std::cout << "msg3: " << msg3 << std::endl;
    std::cout << std::endl;

    Plaintext ptx1 = cc->MakeCKKSPackedPlaintext(msg1);
    Plaintext ptx2 = cc->MakeCKKSPackedPlaintext(msg2);
    Plaintext ptx3 = cc->MakeCKKSPackedPlaintext(msg3);
    Plaintext ptx4 = cc->MakeCKKSPackedPlaintext(msk1);
    Plaintext ptx5 = cc->MakeCKKSPackedPlaintext(msk2);
    Plaintext ptx6 = cc->MakeCKKSPackedPlaintext(msk3);

    auto ctx1 = cc->Encrypt(keyPair.publicKey, ptx1);
    auto ctx2 = cc->Encrypt(keyPair.publicKey, ptx2);
    auto ctx3 = cc->Encrypt(keyPair.publicKey, ptx3);

    ctx1 = cc->EvalLazyRotate(ctx1, 1); 
    ctx2 = cc->EvalLazyRotate(ctx2, 2); 
    ctx3 = cc->EvalLazyRotate(ctx3, 3); 
    auto ctx_add = cc->EvalLazyAdd(ctx1, ctx2);
    ctx_add = cc->EvalLazyAdd(ctx_add, ctx3);
    ctx_add = cc->EvalMult(ctx_add, ptx4);
    
    ctx1 = cc->EvalBatchedKS(ctx1);
    ctx2 = cc->EvalBatchedKS(ctx2);
    ctx3 = cc->EvalBatchedKS(ctx3);
    ctx_add = cc->EvalBatchedKS(ctx_add);
    
    Plaintext decrypted_ptx1, decrypted_ptx2, decrypted_ptx3, decrypted_ptx4;

    cc->Decrypt(keyPair.secretKey, ctx1, &decrypted_ptx1);
    decrypted_ptx1->SetLength(4);
    std::vector<double> decrypted_msg1 = decrypted_ptx1->GetRealPackedValue();
    std::cout << "msg1: " << decrypted_msg1 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx2, &decrypted_ptx2);
    decrypted_ptx2->SetLength(4);
    std::vector<double> decrypted_msg2 = decrypted_ptx2->GetRealPackedValue();
    std::cout << "msg2: " << decrypted_msg2 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx3, &decrypted_ptx3);
    decrypted_ptx3->SetLength(4);
    std::vector<double> decrypted_msg3 = decrypted_ptx3->GetRealPackedValue();
    std::cout << "msg3: " << decrypted_msg3 << std::endl;
    
    cc->Decrypt(keyPair.secretKey, ctx_add, &decrypted_ptx4);
    decrypted_ptx4->SetLength(4);
    std::vector<double> decrypted_msg4 = decrypted_ptx4->GetRealPackedValue();
    std::cout << "add: " << decrypted_msg4 << std::endl;

    std::cout << std::endl;
}
void BFVTest(){
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(2);
    parameters.SetRingDim(1<<7);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetKeySwitchTechnique(BATCHED); // new key switching technique
    
    std::cout << "Using Key Switching Technique: " << parameters.GetKeySwitchTechnique() << std::endl;
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    std::cout << "BFV scheme is using ring dimension "
              << cc->GetRingDimension() << std::endl
              << std::endl;
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);
    cc->Enable(KEYSWITCH);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    // cc->EvalLazyRotateKeyGen(keyPair.secretKey, {1,2,3});
    cc->EvalRotateKeyGen(keyPair.secretKey, {1,2,3});

    std::vector<int64_t> msg1 = {1, 2, 3, 4};
    std::vector<int64_t> msg2 = {5, 6, 7, 8};
    std::vector<int64_t> msg3 = {9, 10, 11, 12};

    std::vector<int64_t> msk1 = {1, 0, 0, 0};
    std::vector<int64_t> msk2 = {0, 1, 0, 0};
    std::vector<int64_t> msk3 = {0, 0, 1, 0};

    std::cout << "msg1: " << msg1 << std::endl;
    std::cout << "msg2: " << msg2 << std::endl;
    std::cout << "msg3: " << msg3 << std::endl;
    std::cout << std::endl;

    Plaintext ptx1 = cc->MakePackedPlaintext(msg1);
    Plaintext ptx2 = cc->MakePackedPlaintext(msg2);
    Plaintext ptx3 = cc->MakePackedPlaintext(msg3);
    Plaintext ptx4 = cc->MakePackedPlaintext(msk1);
    Plaintext ptx5 = cc->MakePackedPlaintext(msk2);
    Plaintext ptx6 = cc->MakePackedPlaintext(msk3);

    std::cout << "encoding" << std::endl;
    
    auto ctx1 = cc->Encrypt(keyPair.publicKey, ptx1);
    auto ctx2 = cc->Encrypt(keyPair.publicKey, ptx2);
    auto ctx3 = cc->Encrypt(keyPair.publicKey, ptx3);
    
    std::cout << "encrypt" << std::endl;
    
    ctx1 = cc->EvalRotate(ctx1, 1); 
    ctx2 = cc->EvalRotate(ctx2, 2); 
    ctx3 = cc->EvalRotate(ctx3, 3); 
    // ctx1 = cc->EvalLazyRotate(ctx1, 1); 
    // ctx2 = cc->EvalLazyRotate(ctx2, 2); 
    // ctx3 = cc->EvalLazyRotate(ctx3, 3); 
    
    std::cout << "lazy rotate" << std::endl;
    
    // auto ctx_add = cc->EvalLazyAdd(ctx1, ctx2);
    // ctx_add = cc->EvalLazyAdd(ctx_add, ctx3);
    
    // std::cout << "lazy add" << std::endl;
    // ctx_add = cc->EvalMult(ctx_add, ptx4);
    // std::cout << "lazy mult" << std::endl;
    
    // ctx1 = cc->EvalBatchedKS(ctx1);
    // ctx2 = cc->EvalBatchedKS(ctx2);
    // ctx3 = cc->EvalBatchedKS(ctx3);
    // ctx_add = cc->EvalBatchedKS(ctx_add);
    
    // std::cout << "ks" << std::endl;
    Plaintext decrypted_ptx1, decrypted_ptx2, decrypted_ptx3, decrypted_ptx4;

    cc->Decrypt(keyPair.secretKey, ctx1, &decrypted_ptx1);
    decrypted_ptx1->SetLength(4);
    std::vector<int64_t> decrypted_msg1 = decrypted_ptx1->GetPackedValue();
    std::cout << "msg1: " << decrypted_msg1 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx2, &decrypted_ptx2);
    decrypted_ptx2->SetLength(4);
    std::vector<int64_t> decrypted_msg2 = decrypted_ptx2->GetPackedValue();
    std::cout << "msg2: " << decrypted_msg2 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx3, &decrypted_ptx3);
    decrypted_ptx3->SetLength(4);
    std::vector<int64_t> decrypted_msg3 = decrypted_ptx3->GetPackedValue();
    std::cout << "msg3: " << decrypted_msg3 << std::endl;
    
    // cc->Decrypt(keyPair.secretKey, ctx_add, &decrypted_ptx4);
    // decrypted_ptx4->SetLength(4);
    // std::vector<int64_t> decrypted_msg4 = decrypted_ptx4->GetPackedValue();
    // std::cout << "add: " << decrypted_msg4 << std::endl;

    std::cout << std::endl;
}
void BGVTest(){
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(2);
    parameters.SetRingDim(1<<7);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetKeySwitchTechnique(BATCHED); // new key switching technique
    
    std::cout << "Using Key Switching Technique: " << parameters.GetKeySwitchTechnique() << std::endl;
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    std::cout << "BGV scheme is using ring dimension "
              << cc->GetRingDimension() << std::endl
              << std::endl;
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);
    cc->Enable(KEYSWITCH);

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    // cc->EvalLazyRotateKeyGen(keyPair.secretKey, {1,2,3});
    cc->EvalRotateKeyGen(keyPair.secretKey, {1,2,3});

    std::vector<int64_t> msg1 = {1, 2, 3, 4};
    std::vector<int64_t> msg2 = {5, 6, 7, 8};
    std::vector<int64_t> msg3 = {9, 10, 11, 12};

    std::vector<int64_t> msk1 = {1, 0, 0, 0};
    std::vector<int64_t> msk2 = {0, 1, 0, 0};
    std::vector<int64_t> msk3 = {0, 0, 1, 0};

    std::cout << "msg1: " << msg1 << std::endl;
    std::cout << "msg2: " << msg2 << std::endl;
    std::cout << "msg3: " << msg3 << std::endl;
    std::cout << std::endl;

    Plaintext ptx1 = cc->MakePackedPlaintext(msg1);
    Plaintext ptx2 = cc->MakePackedPlaintext(msg2);
    Plaintext ptx3 = cc->MakePackedPlaintext(msg3);
    Plaintext ptx4 = cc->MakePackedPlaintext(msk1);
    Plaintext ptx5 = cc->MakePackedPlaintext(msk2);
    Plaintext ptx6 = cc->MakePackedPlaintext(msk3);

    std::cout << "encoding" << std::endl;
    
    auto ctx1 = cc->Encrypt(keyPair.publicKey, ptx1);
    auto ctx2 = cc->Encrypt(keyPair.publicKey, ptx2);
    auto ctx3 = cc->Encrypt(keyPair.publicKey, ptx3);
    
    std::cout << "encrypt" << std::endl;
    
    ctx1 = cc->EvalRotate(ctx1, 1); 
    ctx2 = cc->EvalRotate(ctx2, 2); 
    ctx3 = cc->EvalRotate(ctx3, 3); 
    // ctx1 = cc->EvalLazyRotate(ctx1, 1); 
    // ctx2 = cc->EvalLazyRotate(ctx2, 2); 
    // ctx3 = cc->EvalLazyRotate(ctx3, 3); 
    
    std::cout << "lazy rotate" << std::endl;
    
    // auto ctx_add = cc->EvalLazyAdd(ctx1, ctx2);
    // ctx_add = cc->EvalLazyAdd(ctx_add, ctx3);
    
    // std::cout << "lazy add" << std::endl;
    // ctx_add = cc->EvalMult(ctx_add, ptx4);
    // std::cout << "lazy mult" << std::endl;
    
    // ctx1 = cc->EvalBatchedKS(ctx1);
    // ctx2 = cc->EvalBatchedKS(ctx2);
    // ctx3 = cc->EvalBatchedKS(ctx3);
    // ctx_add = cc->EvalBatchedKS(ctx_add);
    
    // std::cout << "ks" << std::endl;
    Plaintext decrypted_ptx1, decrypted_ptx2, decrypted_ptx3, decrypted_ptx4;

    cc->Decrypt(keyPair.secretKey, ctx1, &decrypted_ptx1);
    decrypted_ptx1->SetLength(4);
    std::vector<int64_t> decrypted_msg1 = decrypted_ptx1->GetPackedValue();
    std::cout << "msg1: " << decrypted_msg1 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx2, &decrypted_ptx2);
    decrypted_ptx2->SetLength(4);
    std::vector<int64_t> decrypted_msg2 = decrypted_ptx2->GetPackedValue();
    std::cout << "msg2: " << decrypted_msg2 << std::endl;
  
    cc->Decrypt(keyPair.secretKey, ctx3, &decrypted_ptx3);
    decrypted_ptx3->SetLength(4);
    std::vector<int64_t> decrypted_msg3 = decrypted_ptx3->GetPackedValue();
    std::cout << "msg3: " << decrypted_msg3 << std::endl;
    
    // cc->Decrypt(keyPair.secretKey, ctx_add, &decrypted_ptx4);
    // decrypted_ptx4->SetLength(4);
    // std::vector<int64_t> decrypted_msg4 = decrypted_ptx4->GetPackedValue();
    // std::cout << "add: " << decrypted_msg4 << std::endl;

    std::cout << std::endl;
    
}

int main() {
    CKKSTest();
    // BFVTest();
    // BGVTest();
    return 0;
}
