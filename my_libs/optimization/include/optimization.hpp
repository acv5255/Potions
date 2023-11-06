#pragma once

#include <stdexcept>
#include <functional>
#include <optional>
#include "armadillo"

using std::function;
using std::optional;
using arma::Col;
using arma::Mat;

const unsigned int MAX_ITERATIONS = 50;

template<typename T>
Mat<T> jacobian(const function<Col<T>(Col<T>)>& func, const Col<T>& x);

template<typename T>
optional<Col<T>> root(const function<Col<T>(Col<T>)>& func, const function<Mat<T>(Col<T>)>& jac, const Col<T>& x0);
