#include "utils-matrices.h"
#include "utils-ptxt.h"
#include "bench_stages.h"
#include "rotation_collector_lazy.h"
#include <cassert>

using namespace benchstages;

static Ciphertext<DCRTPoly> makeZeroCT(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    int slots)
{
    std::vector<double> zero(static_cast<size_t>(slots), 0.0);
    auto pt = cc->MakeCKKSPackedPlaintext(zero);
    return cc->Encrypt(pk, pt);
}
static inline void observeIfBoundary2(RotationKeyCollectorLazy& rk,
    const Ciphertext<DCRTPoly>& ct,
    int& rot_count,
    bool last)
    {
        ++rot_count;
        if (rot_count == 2 || last) {
            rk.observeAutoIndices(ct->GetElementKeyIndexVector());
            rot_count = 0;
        }
    }
static inline void observeIfBoundary3(RotationKeyCollectorLazy& rk,
    const Ciphertext<DCRTPoly>& ct,
    int& rot_count,
    bool last)
    {
        ++rot_count;
        if (rot_count == 3 || last) {
            rk.observeAutoIndices(ct->GetElementKeyIndexVector());
            rot_count = 0;
        }
    }
    
    Ciphertext<DCRTPoly> maskRow(
        Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t rowIndex
)
{
    assert(rowIndex >= 0 && rowIndex < matrixSize && "Invalid row index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * rowIndex + i] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> maskColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    const size_t columnIndex
)
{
    assert(columnIndex >= 0 && columnIndex < matrixSize && "Invalid column index");

    std::vector<double> mask(matrixSize * matrixSize, 0.0);
    for (size_t i = 0; i < matrixSize; i++)
        mask[matrixSize * i + columnIndex] = 1.0;
    
    return c * mask;
}


Ciphertext<DCRTPoly> replicateRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    return c;
}

// 3 LazyRot -> BatchedKS
Ciphertext<DCRTPoly> replicateRowLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();

    auto tmp1 = cc->EvalLazyRotate(c, - (1 << (LOG2(matrixSize))));
    auto tmp2 = cc->EvalLazyAdd(c, tmp1);
    rot_count++;
    for (size_t i = 1; i < LOG2(matrixSize); i++){
        tmp1 = cc->EvalLazyRotate(tmp2, - (1 << (LOG2(matrixSize) + i)));
        cc->EvalLazyAddInPlace(tmp2, tmp1);
        rot_count++;

        const bool last = (i ==  LOG2(matrixSize) - 1);
    if (rot_count == 2 || last) {
            tmp2 = cc->EvalBatchedKS(tmp2);
            rot_count = 0;
        }
    }
    return tmp2;
}
void replicateRowLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;
    for (size_t i = 0; i < LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, - (1 << (LOG2(matrixSize) + i)));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) - 1);
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}


Ciphertext<DCRTPoly> replicateColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << i);

    return c;
}

Ciphertext<DCRTPoly> replicateColumnLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();
    for (size_t i = 0; i < LOG2(matrixSize); i++){
        auto tmp = cc->EvalLazyRotate(c, - (1 << i));
        cc->EvalLazyAddInPlace(c, tmp);
        rot_count++;

        const bool last = (i ==  LOG2(matrixSize) - 1);
        if (rot_count == 2 || last) {
            c = cc->EvalBatchedKS(c);
            rot_count = 0;
        }
    }
    return c;
}

void replicateColumnLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;
    for (size_t i = 0; i < LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, - (1 << i));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) - 1);
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}

Ciphertext<DCRTPoly> sumRows(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput,
    const size_t outputRow
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c >> (1 << (LOG2(matrixSize) + i));

    if (maskOutput)
        c = maskRow(c, matrixSize, outputRow);

    return c;
}
Ciphertext<DCRTPoly> sumRowsLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput,
    const size_t outputRow
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();
    for (size_t i = 0; i < LOG2(matrixSize); i++){
        auto tmp = cc->EvalLazyRotate(c, - (1 << (LOG2(matrixSize) + i)));
        cc->EvalLazyAddInPlace(c, tmp);
        rot_count++;

        const bool last = (i ==  LOG2(matrixSize) - 1);
        if (rot_count == 2 || last) {
            c = cc->EvalBatchedKS(c);
            rot_count = 0;
        }
    }
    if (maskOutput)
        c = maskRow(c, matrixSize, outputRow);

    return c;
}
void sumRowsLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk,
    bool maskOutput,
    const size_t outputRow)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;
    for (size_t i = 0; i < LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, - (1 << (LOG2(matrixSize) + i)));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) - 1);
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}


Ciphertext<DCRTPoly> sumColumns(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 0; i < LOG2(matrixSize); i++)
        c += c << (1 << i);
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);
    
    return c;
}
Ciphertext<DCRTPoly> sumColumnsLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();
    for (size_t i = 0; i < LOG2(matrixSize); i++){
        auto tmp = cc->EvalLazyRotate(c, (1 << i));
        cc->EvalLazyAddInPlace(c, tmp);
        rot_count++;

        const bool last = (i ==  LOG2(matrixSize) - 1);
        if (rot_count == 2 || last) {
            c = cc->EvalBatchedKS(c);
            rot_count = 0;
        }
    }
    
    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);
    
    return c;
}
void sumColumnsLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk,
    bool maskOutput)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;
    for (size_t i = 0; i < LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, (1 << i));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) - 1);
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}

Ciphertext<DCRTPoly> transposeRow(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c >> (matrixSize * (matrixSize - 1) / (1 << i));

    if (maskOutput)
        c = maskColumn(c, matrixSize, 0);

    return c;
}

Ciphertext<DCRTPoly> transposeRowLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();

    auto tmp1 = cc->EvalLazyRotate(c, - (matrixSize * (matrixSize - 1) / (1 << 1)));
    auto tmp2 = cc->EvalLazyAdd(c, tmp1);
    rot_count++;
    for (size_t i = 2; i <= LOG2(matrixSize); i++){
        tmp1 = cc->EvalLazyRotate(tmp2, - (matrixSize * (matrixSize - 1) / (1 << i)));
        cc->EvalLazyAddInPlace(tmp2, tmp1);
      
        rot_count++;
        const bool last = (i ==  LOG2(matrixSize));
        if (rot_count == 2 || last) {
            tmp2 = cc->EvalBatchedKS(tmp2);
            rot_count = 0;
        }
    }
    
    if (maskOutput)
        tmp2 = maskColumn(tmp2, matrixSize, 0);

    return tmp2;
}
void transposeRowLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk,
    bool maskOutput)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;

    for (size_t i = 1; i <= LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, - (matrixSize * (matrixSize - 1) / (1 << i)));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) );
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}


Ciphertext<DCRTPoly> transposeColumn(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    for (size_t i = 1; i <= LOG2(matrixSize); i++)
        c += c << (matrixSize * (matrixSize - 1) / (1 << i));

    if (maskOutput)
        c = maskRow(c, matrixSize, 0);

    return c;
}

Ciphertext<DCRTPoly> transposeColumnLazy(
    Ciphertext<DCRTPoly> c,
    const size_t matrixSize,
    bool maskOutput
)
{
    int rot_count = 0;
    CryptoContext<DCRTPoly> cc = c->GetCryptoContext();
    for (size_t i = 0; i < LOG2(matrixSize); i++){
        auto tmp = cc->EvalLazyRotate(c, (matrixSize * (matrixSize - 1) / (1 << i)));
        cc->EvalLazyAddInPlace(c, tmp);
        rot_count++;

        const bool last = (i ==  LOG2(matrixSize) - 1);
        if (rot_count == 2 || last) {
            c = cc->EvalBatchedKS(c);
            rot_count = 0;
        }
    }

    if (maskOutput)
        c = maskRow(c, matrixSize, 0);

    return c;
}
void transposeColumnLazyPlan(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    const size_t matrixSize,
    RotationKeyCollectorLazy& rk,
    bool maskOutput)
{
    const int slots = static_cast<int>(matrixSize * matrixSize);
    auto ct = makeZeroCT(cc, pk, slots);
    int rot_count = 0;
    for (size_t i = 0; i < LOG2(matrixSize); i++) {
        auto tmp = cc->EvalLazyRotate(ct, (matrixSize * (matrixSize - 1) / (1 << i)));
        cc->EvalLazyAddInPlace(ct, tmp);
        const bool last = (i == LOG2(matrixSize) - 1);
        observeIfBoundary2(rk, ct, rot_count, last);
    }
}

std::vector<int32_t> getRotationIndices
(
    const size_t matrixSize
)
{
    std::vector<int32_t> indices;

    int32_t index;
    for (size_t i = 0; i < LOG2(matrixSize); i++)
    {
        index = 1 << i;
        indices.push_back(index);   // sumColumns
        indices.push_back(-index);  // replicateColumn

        index = 1 << (LOG2(matrixSize) + i);
        indices.push_back(-index);  // replicateRow, sumRows

        index = (matrixSize * (matrixSize - 1) / (1 << (i + 1)));
        indices.push_back(index);   // transposeColumn
        indices.push_back(-index);  // transposeRow
    }

    return indices;
}


// int main()
// {

//     const usint matrixSize              = 8;

//     const usint integralPrecision       = 10;
//     const usint decimalPrecision        = 50;
//     const usint multiplicativeDepth     = 2;
//     const usint numSlots                = matrixSize * matrixSize;
//     const bool enableBootstrap          = false;
//     const usint ringDim                 = 0;
//     const bool verbose                  = true;

//     std::vector<int32_t> indices = getRotationIndices(matrixSize);

//     CryptoContext<DCRTPoly> cryptoContext = generateCryptoContext(
//         integralPrecision,
//         decimalPrecision,
//         multiplicativeDepth,
//         numSlots,
//         enableBootstrap,
//         ringDim,
//         verbose
//     );

//     KeyPair<DCRTPoly> keyPair = keyGeneration(
//         cryptoContext,
//         indices,
//         numSlots,
//         enableBootstrap,
//         verbose
//     );

//     std::vector<double> matrix(matrixSize * matrixSize);
//     for (size_t i = 0; i < matrixSize; i++)
//         for (size_t j = 0; j < matrixSize; j++)
//             matrix[i * matrixSize + j] = i * matrixSize + j;

//     std::cout << "Matrix: " << vector2matrix(matrix, matrixSize) << std::endl;

//     Ciphertext<DCRTPoly> matrixC = cryptoContext->Encrypt(
//         keyPair.publicKey,
//         cryptoContext->MakeCKKSPackedPlaintext(matrix)
//     );

//     auto start = std::chrono::high_resolution_clock::now();

//     Ciphertext<DCRTPoly> resultC = matrixC;
//     // Ciphertext<DCRTPoly> resultC = maskRow(matrixC, matrixSize, 0);
//     // resultC = replicateRow(resultC, matrixSize);
//     // resultC = transposeRow(resultC, matrixSize, true);
//     // Ciphertext<DCRTPoly> resultC = maskColumn(matrixC, matrixSize, 0);
//     // resultC = replicateColumn(resultC, matrixSize);
//     // resultC = transposeColumn(resultC, matrixSize, true);
//     // Ciphertext<DCRTPoly> resultC = sumRows(matrixC, matrixSize, index);
//     // Ciphertext<DCRTPoly> resultC = sumColumns(matrixC, matrixSize);

//     auto end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = end - start;
//     std::cout << elapsed_seconds.count() << "s" << std::endl;

//     Plaintext resultP;
//     cryptoContext->Decrypt(keyPair.secretKey, resultC, &resultP);
//     resultP->SetLength(matrixSize * matrixSize);

//     std::vector<std::vector<double>> result = vector2matrix(
//         resultP->GetRealPackedValue(),
//         matrixSize
//     );
//     std::cout << "Result: " << result << std::endl;

//     return 0;

// }
