/*
 * This code implements algorithms from:
 * "Efficient Ranking, Order Statistics, and Sorting under CKKS"
 * by Federico Mazzone, Maarten H. Everts, Florian Hahn, and Andreas Peter
 * (https://doi.org/10.48550/arXiv.2412.15126)
 *
 * this implementation is from:
 * https://github.com/FedericoMazzone/openfhe-statistics
 * Copyright (c) 2024 Federico Mazzone
 * Licensed under BSD 2-Clause License
 */
#pragma once

#include "openfhe.h"
#include "rotation_collector_lazy.h"
using namespace lbcrypto;

/**
 * @brief Sorts the elements of a vector in ascending order.
 * 
 * This function sorts the elements of the input vector `vec` in ascending
 * order.
 * 
 * @param vec The input vector of doubles to be sorted.
 * @return std::vector<double> The sorted form of `vec`.
 */
std::vector<double> sort(
    std::vector<double> vec
);


/**
 * @brief Sorts a ciphertext vector.
 * 
 * This function sorts a ciphertext vector `c`.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param degreeI The degree of the indicator function's approximation.
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sort(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    double leftBoundC,
    double rightBoundC,
    uint32_t degreeC,
    uint32_t degreeI
);

/**
 * @brief Sorts a ciphertext vector.
 * 
 * This function sorts a ciphertext vector `c`. It uses the fg approximation of
 * the sign function.
 * 
 * @param c The ciphertext vector to sort.
 * @param vectorLength The length of the vector.
 * @param dg_c The composition degree of g for cmp (reduce the input gap).
 * @param df_c The composition degree of f for cmp (reduce the output error).
 * @param dg_i The composition degree of g for ind (reduce the input gap).
 * @param df_i The composition degree of f for ind (reduce the output error).
 * @return Ciphertext<DCRTPoly> A ciphertext representing the sorted vector.
 */
Ciphertext<DCRTPoly> sortFG(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
);

Ciphertext<DCRTPoly> sortFGLazy(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i
);
Ciphertext<DCRTPoly> sortFGLazyDebug(
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i,
    const PrivateKey<DCRTPoly>& sk,
    size_t preview = 4,          // 출력 미리보기 개수
    int precision = 3            // 소수점 자리
);


void Plan_sortFGLazy(
    const CryptoContext<DCRTPoly>& cc,
    const PublicKey<DCRTPoly>& pk,
    size_t vectorLength,
    uint32_t dg_c,
    uint32_t df_c,
    uint32_t dg_i,
    uint32_t df_i,
    RotationKeyCollectorLazy& rk
);
