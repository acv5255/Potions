#include "catch2/catch_test_macros.hpp"
#include "optimization.hpp"

using arma::vec;
using arma::mat;

// TEST_CASE("optimization::jacobian") {
//     REQUIRE(false);
// }

TEST_CASE("optimization::root") {
    /*
        Test a few functions:
        0: y = x(x + 25) (root @ 0)
        1: y = e^x - 1         (root @ 0)
        2: y = x(x + 2)(x + 3)(x - 10) (root @ 0)
     */
    auto f = [&] (const vec& x) -> vec {
        vec y = arma::zeros(x.size());

        y[0] = x[0] * (x[0] + 25.0);
        y[1] = std::exp(x[1]) - 1.0;
        y[2] = x[2] * (x[2] + 2.0) * (x[2] + 3) * (x[2] - 10);
        return y;
    };

    vec x0 = arma::ones(3);
    const vec roots_actual = arma::zeros(3);

    auto fprime = [&] (const vec& x) -> mat {
        return jacobian<double>(f, x);
    };

    optional<vec> roots_opt = root<double>(f, fprime, x0);

    REQUIRE(roots_opt.has_value());
    vec roots = roots_opt.value();

    const vec err = roots_actual - roots;
    const double max_err = arma::abs(err).max();

    REQUIRE(max_err < 1e-8);
}