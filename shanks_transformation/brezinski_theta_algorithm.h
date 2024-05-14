﻿/**
 * @file theta_brezinski_algorithm.h
 * @brief This file contains the declaration of the Theta Brezinski Algorithm class.
 */

#pragma once

#include "series_acceleration.h" // Include the series header
//#include <vector> // Include the vector library
#define epsilon 1.0E-160 //precision of calculations. checks numeders to being close to zero

template <typename T, typename K, typename series_templ>
class theta_brezinski_algorithm : public series_acceleration<T, K, series_templ>
{
public:
    /*
     * @brief Parameterized constructor to initialize the Theta Brezinski Algorithm.
     * @authors Yurov P.I. Bezzaborov A.A.
     * @param series The series class object to be accelerated
     */
    theta_brezinski_algorithm(const series_templ& series);

    /**
     * @brief Fast implementation of Theta Brezinski algorithm.
     * Computes the partial sum after the transformation using the Theta Brezinski Algorithm.
     * For more information, see p. 277 10.2-4 in [https://arxiv.org/pdf/math/0306302.pdf]
     * @param n The number of terms in the partial sum.
     * @param order The order of transformation.
     * @return The partial sum after the transformation.
     */
    T operator()(const K n, const int order) const;

protected:
    /**
     * @brief Recursive function to compute theta.
     * Computes the value of theta according to the given parameters.
     * @param n The number of terms in the partial sum.
     * @param k The order of transformation.
     * @return The value of theta.
     */
    T theta(K n, const int k, T S_n, const K j) const;
};

template <typename T, typename K, typename series_templ>
theta_brezinski_algorithm<T, K, series_templ>::theta_brezinski_algorithm(const series_templ& series)
    : series_acceleration<T, K, series_templ>(series) {}

template <typename T, typename K, typename series_templ>
T theta_brezinski_algorithm<T, K, series_templ>::operator()(const K n, const int order) const
{
    if (order % 2 != 0 || order < 0) 
        throw std::domain_error("order should be even number");
    else if (n < 0)
        throw std::domain_error("negative integer in the input");
    else if (n == 0 || order == 0)
        return this->series->S_n(n);
    return theta(n, order, this->series->S_n(n),0);
}


template <typename T, typename K, typename series_templ>
T theta_brezinski_algorithm<T, K, series_templ>::theta(K n, const int k, T S_n, const K j) const
{
    if (k == 1)
    {
        T res = 1/this->series->operator()(n+j+1);
        if (!std::isfinite(res)) throw std::overflow_error("division by zero");
        return res;
    }

    for (K tmp = n+1; tmp <= n + j; tmp++) {
        S_n += this->series->operator()(tmp);
    }
    n += j;

    if (k == 0)
    {
        return S_n;
    }
    else if (k % 2 == 1)
    {
        const T delta = theta(n, k - 1, S_n, 0) - theta(n, k - 1, S_n, 1); // Δυ_2k^(n)

        if(delta > -epsilon && delta < epsilon) throw std::overflow_error("division by zero");

        return theta(n, k - 2,S_n ,1) + T(1) / delta; // υ_(2k+1)^(n)=υ_(2k-1)^(n+1) + 1/(Δυ_2k^(n)
        
    }
    else
    { // k % 2 == 0
        const T delta2 = theta(n, k - 1, S_n, 0) - 2 * theta(n, k - 1,S_n,1) + theta(n, k - 1,S_n,2); // Δ^2 υ_(2k+1)^(n)

        if (delta2 > -epsilon && delta2 < epsilon) throw std::overflow_error("division by zero");

        const T delta_n = theta(n, k - 2, S_n, 1) - theta(n , k - 2, S_n,2); // Δυ_2k^(n+1) 
        const T delta_n1 = theta(n, k - 1, S_n,1) - theta(n, k - 1,S_n,2); // Δυ_(2k+1)^(n+1)
        
        return theta(n, k - 2, S_n, 1) + (delta_n * delta_n1) / delta2; // υ_(2k+2)^(n)=υ_2k^(n+1)+((Δυ_2k^(n+1))*(Δυ_(2k+1)^(n+1)))/(Δ^2 υ_(2k+1)^(n)
    }
}