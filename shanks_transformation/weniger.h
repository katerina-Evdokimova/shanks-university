/**
* @file weniger.h
* @brief This files contains the definition of Weniger-transformation
*/
#pragma once
#define DEF_UNDEFINED_SUM 0

#include "series_acceleration.h" // Include the series header
#include <vector>
#include <iostream>

template<typename T, typename K, typename series_templ>
class weniger_algorithm : public series_acceleration<T, K, series_templ>
{
protected:

	/**
	* @brief Default function to calculate W-tranformation.
	* For more information see "Joint use of the Weniger transformation and hyperasymptotics for accurate asymptotic evaluations of a class of saddle-point integrals"
	* @authors Yurov P.I. Bezzaborov A.A.
	* @param k The number of terms in the partial sum.
	* @param n the order of transformation
	* @param remainder_func functor, whose returning w_n for t,u or v transformation
	* @return The partial sum after the transformation.
	*/


public:

	
	/**
	* @brief Weniger class template for derivations
	* @authors Yurov P.I. Bezzaborov A.A.
	* @tparam T The type of the elements in the series, K The type of enumerating integer, series_templ is the type of series whose convergence we accelerate
	*/

	weniger_algorithm(const series_templ& series);

	/**
	* @brief Abstract method for the inherited classes below
	* Computes the partial sum after the transformation using the epsilon Weniger Algorithm.
	* @param k The number of terms in the partial sum.
	* @param n The order of transformation.
	* @return The partial sum after the transformation.
	*/

	T operator()(const K k, const int n) const;

};

template<typename T, typename K, typename series_templ>
weniger_algorithm<T, K, series_templ>::weniger_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template<typename T, typename K, typename series_templ>
T weniger_algorithm<T, K, series_templ>::operator()(const K k, const int n) const {

	T numerator = T(0), denominator = T(0);
	T a_n, rest;
	T coef = T(1);


	T binomial_coef = weniger_algorithm<T, K, series_templ>::series->binomial_coefficient(k, 0);
	T S_n = weniger_algorithm<T, K, series_templ>::series->S_n(0);

	for (int m = 0; m < k - 1; ++m) {
		coef *= (1 + m);

	}

	for (int j = 0; j <= k; ++j) {

		rest = weniger_algorithm<T, K, series_templ>::series->minus_one_raised_to_power_n(j) * binomial_coef;
		binomial_coef = binomial_coef * (k - j) / (j + 1);

		rest = rest * coef;

		coef = coef / (1 + j) * (1 + j + 1 + k - 2);

		a_n = 1 / weniger_algorithm<T, K, series_templ>::series->operator()(j + 1);

		numerator += rest * S_n * a_n;

		S_n += weniger_algorithm<T, K, series_templ>::series->operator()(j + 1);

		denominator += rest * a_n;

	}

	numerator /= denominator;

	if (!std::isfinite(numerator))
		throw std::overflow_error("division by zero");

	return numerator;
}