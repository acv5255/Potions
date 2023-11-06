#pragma once

#include <stdexcept>
#include <functional>
#include <optional>
#include "armadillo"

using std::function;
using std::optional;
using arma::Col;
using arma::Mat;

const double TOLERANCE = 1e-10;
const unsigned int MAX_ITERATIONS = 50;

template<typename T>
Mat<T> jacobian(const function<Col<T>(Col<T>)>& func, const Col<T>& x);

template<typename T>
optional<Col<T>> root(const function<Col<T>(Col<T>)>& func, const function<Mat<T>(Col<T>)>& jac, const Col<T>& x0);

class NotImplemented : public std::logic_error {
    public:
        NotImplemented() : std::logic_error("Function not yet implemented") { };
};