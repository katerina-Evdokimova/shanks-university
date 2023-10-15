#pragma once
#include <functional>
#include <iostream>
#include <exception>
#include <math.h>

template <typename T>
class series_acceleration
{
public:
	series_acceleration();
	series_acceleration(std::function<T(T, int)> series, T x);
	virtual ~series_acceleration() = 0;
	void print_s_n(int n);
	void print_t_n(int n);
protected:
	virtual T transform(int n);
	std::function<T(T, int)> series;
	T x;
	T S_n(int n);
};

template <typename T>
series_acceleration<T>::series_acceleration()
{
	series = [](T x, int n) -> T {return 0; };
}

template <typename T>
series_acceleration<T>::series_acceleration(std::function<T(T, int)> series, T x)
{
	this->series = series;
	this->x = x;
}

template <typename T>
series_acceleration<T>::~series_acceleration()
{

}

template <typename T>
void series_acceleration<T>::print_s_n(int n)
{
	T s_n = 0;	//its ineffective to use s_n method if n is big
	for (int i = 0; i <= n; ++i)	s_n += series(x, i);
	std::cout << "S_" << n << " : " << s_n << std::endl;
}

template <typename T>
void series_acceleration<T>::print_t_n(int n)
{
	std::cout << "T_" << n << " : " << transform(n) << std::endl;
}

template <typename T>
T series_acceleration<T>::S_n(int n)
{
	if (n < 0)	throw std::out_of_range("negative n"); //to do - specify exception
	return n == 0 ? series(1,0) : S_n(n - 1) + series(x, n);
}

template <typename T>
T series_acceleration<T>::transform(int n) { return 0; }