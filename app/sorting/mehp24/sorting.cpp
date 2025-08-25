/*
 * This code implements algorithms from:
 * "Efficient Ranking, Order Statistics, and Sorting under CKKS"
 * by Federico Mazzone, Maarten H. Everts, Florian Hahn, and Andreas Peter
 * (https://doi.org/10.48550/arXiv.2412.15126)
 *
 * this implementation are from:
 * https://github.com/FedericoMazzone/openfhe-statistics
 * Copyright (c) 2024 Federico Mazzone
 * Licensed under BSD 2-Clause License
 */

#include "sorting.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include "bench_stages.h"
#include "rotation_collector_lazy.h"
#include <cassert>
#include <omp.h>

using namespace benchstages;

std::vector<double> sort(
    std::vector<double> vec
)
{
    std::sort(vec.begin(), vec.end());
    return vec;
}


Ciphertext<DCRTPoly> sort(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
)
{
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    Ciphertext<DCRTPoly> VC = replicateColumn(transposeRow(c, vectorLength, true), vectorLength);
    Ciphertext<DCRTPoly> C = compare(
        VR, VC,
        leftBoundC, rightBoundC,
        degreeC
    );
    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength);

    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 0.5;
    Ciphertext<DCRTPoly> M = indicator(
        R + subMask,
        -0.5, 0.5,
        -1.0 * vectorLength, 1.0 * vectorLength,
        degreeI
    );

    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);

    return S;
}

Ciphertext<DCRTPoly> sortFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
)
{
    // Stage 1: replicateRow
    benchstages::stages_begin(benchstages::ST_REPLICATE_ROW);
    Ciphertext<DCRTPoly> VR = replicateRow(c, vectorLength);
    benchstages::stages_end(benchstages::ST_REPLICATE_ROW);

    // Stage 2: transposeRow
    benchstages::stages_begin(benchstages::ST_TRANSPOSE_ROW);
    Ciphertext<DCRTPoly> TR = transposeRow(c, vectorLength, true);
    benchstages::stages_end(benchstages::ST_TRANSPOSE_ROW);

    // Stage 3: replicateColumn
    benchstages::stages_begin(benchstages::ST_REPLICATE_COL);
    Ciphertext<DCRTPoly> VC = replicateColumn(TR, vectorLength);
    benchstages::stages_end(benchstages::ST_REPLICATE_COL);

    // Stage 4: compareAdv
    benchstages::stages_begin(benchstages::ST_COMPARE_ADV);
    Ciphertext<DCRTPoly> C = compareAdv(
        VR, VC,
        dg_c, df_c
    );
    benchstages::stages_end(benchstages::ST_COMPARE_ADV);
    // std::cout << "C levels: " << C->GetLevel() << std::endl;

    // Stage 5: sumRows
    benchstages::stages_begin(benchstages::ST_SUM_ROWS);
    Ciphertext<DCRTPoly> R = sumRows(C, vectorLength);
    benchstages::stages_end(benchstages::ST_SUM_ROWS);

    // Stage 6: indicatorAdv
    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 0.5;
    benchstages::stages_begin(benchstages::ST_INDICATOR);
    Ciphertext<DCRTPoly> M = indicatorAdv(
        R + subMask,
        vectorLength,
        dg_i, df_i
    );
    benchstages::stages_end(benchstages::ST_INDICATOR);
    // std::cout << "M levels: " << M->GetLevel() << std::endl;

    // Stage 7: sumColumns
    benchstages::stages_begin(benchstages::ST_SUM_COLS);
    Ciphertext<DCRTPoly> S = sumColumns(M * VR, vectorLength);
    benchstages::stages_end(benchstages::ST_SUM_COLS);
    // std::cout << "S levels: " << S->GetLevel() << std::endl;

    return S;
}

Ciphertext<DCRTPoly> sortFGLazy(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
)
{
    // Stage 1: replicateRowLazy
    benchstages::stages_begin(benchstages::ST_REPLICATE_ROW);
    Ciphertext<DCRTPoly> VR = replicateRowLazy(c, vectorLength);
    auto vr = VR->Clone();
    benchstages::stages_end(benchstages::ST_REPLICATE_ROW);
    
    // Stage 2: transposeRowLazy
    benchstages::stages_begin(benchstages::ST_TRANSPOSE_ROW);
    Ciphertext<DCRTPoly> TR = transposeRowLazy(c, vectorLength, true);
    benchstages::stages_end(benchstages::ST_TRANSPOSE_ROW);
    
    // Stage 3: replicateColumnLazy
    benchstages::stages_begin(benchstages::ST_REPLICATE_COL);
    Ciphertext<DCRTPoly> VC = replicateColumnLazy(TR, vectorLength);
    benchstages::stages_end(benchstages::ST_REPLICATE_COL);
    
    // Stage 4: compareAdv
    benchstages::stages_begin(benchstages::ST_COMPARE_ADV);
    Ciphertext<DCRTPoly> C = compareAdv(
        VR, VC,
        dg_c, df_c
    );
    benchstages::stages_end(benchstages::ST_COMPARE_ADV);
    // std::cout << "C levels: " << C->GetLevel() << std::endl;
    
    // Stage 5: sumRowsLazy
    benchstages::stages_begin(benchstages::ST_SUM_ROWS);
    Ciphertext<DCRTPoly> R = sumRowsLazy(C, vectorLength);
    benchstages::stages_end(benchstages::ST_SUM_ROWS);

    // Stage 6: indicatorAdv
    std::vector<double> subMask(vectorLength * vectorLength);
    for (size_t i = 0; i < vectorLength; i++)
        for (size_t j = 0; j < vectorLength; j++)
            subMask[i * vectorLength + j] = -1.0 * i - 0.5;
    benchstages::stages_begin(benchstages::ST_INDICATOR);
    Ciphertext<DCRTPoly> M = indicatorAdv(
        R + subMask,
        vectorLength,
        dg_i, df_i
    );
    benchstages::stages_end(benchstages::ST_INDICATOR);
    M->SetElementKeyIndexVector({0,1});
    // std::cout << "M levels: " << M->GetLevel() << std::endl;

    // Stage 7: sumColumnsLazy
    benchstages::stages_begin(benchstages::ST_SUM_COLS);

    auto MVR = M * VR;

    MVR->SetElementKeyIndexVector({0,1});
    Ciphertext<DCRTPoly> S = sumColumnsLazy(MVR, vectorLength);
    benchstages::stages_end(benchstages::ST_SUM_COLS);
    // std::cout << "S levels: " << S->GetLevel() << std::endl;

    return S;
}

void Plan_sortFGLazy(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    size_t vectorLength,
    uint32_t /*dg_c*/,
    uint32_t /*df_c*/,
    uint32_t /*dg_i*/,
    uint32_t /*df_i*/,
    RotationKeyCollectorLazy& rk)
{
    rk.begin(static_cast<int>(vectorLength * vectorLength));

    // Stage 1
    replicateRowLazyPlan(cc, pk, vectorLength, rk);
    // Stage 2
    transposeRowLazyPlan(cc, pk, vectorLength, rk, /*maskOutput=*/true);
    // Stage 3
    replicateColumnLazyPlan(cc, pk, vectorLength, rk);
    // Stage 4 (compareAdv): no rotations to collect
    // Stage 5
    sumRowsLazyPlan(cc, pk, vectorLength, rk, /*maskOutput=*/false, /*outputRow=*/0);
    // Stage 6 (indicatorAdv): no rotations to collect
    // Stage 7
    sumColumnsLazyPlan(cc, pk, vectorLength, rk, /*maskOutput=*/false);
}





