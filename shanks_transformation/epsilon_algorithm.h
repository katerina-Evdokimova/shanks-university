/**
 * @file epsilon_algorithm.h
 * @brief This file contains the declaration of the Epsilon Algorithm class.
 */

#pragma once
#define DEF_UNDEFINED_SUM 0

#include "series_acceleration.h" // Include the series header
#include <vector> // Include the vector library


 /**
  * @brief Epsilon Algorithm class template.
  * @authors Pashkov B.B.
  * @tparam T The type of the elements in the series, K The type of enumerating integer, series_templ is the type of series whose convergence we accelerate
  */
template <typename T, typename K, typename series_templ>
class epsilon_algorithm : public series_acceleration<T, K, series_templ>
{
public:
	/**
   * @brief Parameterized constructor to initialize the Epsilon Algorithm.
   * @authors Pashkov B.B.
   * @param series The series class object to be accelerated
   */
	epsilon_algorithm(const series_templ& series);

	/**
   * @brief Shanks multistep epsilon algorithm.
   * Computes the partial sum after the transformation using the Epsilon Algorithm.
   * For more information, see p. 5.3.2 in [https://e-maxx.ru/bookz/files/numerical_recipes.pdf]  
   * @param n The number of terms in the partial sum.
   * @param order The order of transformation.
   * @return The partial sum after the transformation.
   */
	T operator()(const K n, const int order) const;
};

template <typename T, typename K, typename series_templ>
epsilon_algorithm<T, K, series_templ>::epsilon_algorithm(const series_templ& series) : series_acceleration<T, K, series_templ>(series) {}

template <typename T, typename K, typename series_templ>
T epsilon_algorithm<T, K, series_templ>::operator()(const K n, const int order) const
{
	int m = 2 * order;
	if (n < 0)
		throw std::domain_error("negative integer in the input");
	else if (n == 0)
		return DEF_UNDEFINED_SUM;
	else if (order == 0)
		return this->series->S_n(n);

	std::vector<T>* e0 = new std::vector<T>(m + n + 1, 0);
	std::vector<T>* e1 = new std::vector<T>(m + n, 0);
	for (int j = m + n; j >= 0; --j)
	{
		(*e0)[j] = this->series->S_n(j);
	}

	int max_ind = m + n;
	for (int i = 0; i < m; ++i)
	{
		for (int j = n-1; j < max_ind; ++j)
		{
			(*e1)[j] += 1.0 / ((*e0)[j + 1] - (*e0)[j]);
		}
		--max_ind;
		std::swap(e0, e1);
		e1->erase(e1->begin());
	}

	const auto result = (*e0)[n - 1];

	if (!std::isfinite(result))
		throw std::overflow_error("division by zero");

	delete e0; delete e1;

	return result;
}