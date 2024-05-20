/**
* @file levin_sidi_merg.h
* @brief This files contains the definition of Drummonds D-transformation with u,t,d,v remainders
*/
#pragma once
#define DEF_UNDEFINED_SUM 0

#include "series_acceleration.h" // Include the series header
#include "levin_sidi_S_algorithm.h"
#include <vector> // Include the vector library



template<typename T, typename K, typename series_templ>
class drummonds_algorithm : public series_acceleration<T, K, series_templ>
{
protected:

	const transform_base<T, K>* remainder_func;

	bool recursive;

	/**
	* @brief Function to calculate D-tranformation directly by formula. For more information see p. 70 9.5-4 [https://arxiv.org/pdf/math/0306302.pdf]
	* @authors Naumov A.
	* @param k The number of terms in the partial sum.
	* @param n the order of transformation
	* @return The partial sum after the transformation.
	*/

	virtual T calculate(const K& k, const int& n) const {

		if (n < 0)
			throw std::domain_error("negative integer in input");


		T numerator = T(0), denominator = T(0);
		T w_n, rest;

		for (int j = 0; j <= k; ++j) {

			rest = this->series->minus_one_raised_to_power_n(j) * this->series->binomial_coefficient(k, j);

			w_n = remainder_func->operator()(n, j, this->series, 1);

			numerator += rest * this->series->S_n(n + j) * w_n;
			denominator += rest * w_n;

		}

		numerator /= denominator;

		if (!std::isfinite(numerator))
			throw std::overflow_error("division by zero");

		return numerator;
	}

	/**
	* @brief Function to calculate D-tranformation using reccurence formula.  For more information see p. 70 9.5-5 [https://arxiv.org/pdf/math/0306302.pdf]
	* @authors Naumov A.
	* @param k The number of terms in the partial sum.
	* @param n the order of transformation
	* @return The partial sum after the transformation.
	*/

	T calculate_rec(const K& k, const int& n) const {

		if (n < 0)
			throw std::domain_error("negative integer in input");

		std::vector<T>* N = new std::vector<T>(k + 1, 0);
		std::vector<T>* D = new std::vector<T>(k + 1, 0);

		for (int i = 0; i < k + 1; i++) {
			(*D)[i] = remainder_func->operator()(0, n + i, this->series);
			(*N)[i] = this->series->S_n(n + i) * (*D)[i];
		}

		for (int i = 1; i <= k; ++i)
			for (int j = 0; j <= k - i; ++j) {

				(*D)[j] = (*D)[j + 1] - (*D)[j];
				(*N)[j] = (*N)[j + 1] - (*N)[j];
			}

		T numerator = (*N)[0] / (*D)[0];

		delete N, D;

		if (!std::isfinite(numerator))
			throw std::overflow_error("division by zero");

		return numerator;
	}

public:

	/**
	* @brief Parameterized constructor to initialize the Drummonds Algorithm.
	* @param series The series class object to be accelerated
	* @param func Remainder function
	* @param recursive How to calculate
	*/

	drummonds_algorithm(const series_templ& series, const transform_base<T, K>* func, bool recursive = false) : series_acceleration<T, K, series_templ>(series), remainder_func(func), recursive(recursive) {}

	~drummonds_algorithm() { delete remainder_func; }

	/**
	* @brief Abstract method for the inherited classes below
	* Computes the partial sum after the transformation using the Drummonds Algorithm.
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @return The partial sum after the transformation.
	*/

	T operator()(const K k, const int n) const {
		if (recursive) return calculate_rec(k, n);
		return calculate(k, n);
	}

};
