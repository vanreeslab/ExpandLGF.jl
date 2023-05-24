#include "gtest/gtest.h"
#include "test.hpp"
#include <array>
#include <algorithm>

using LGFOneUnbounded::real_t;
using std::array;
using std::string;
using std::to_string;
using std::vector;

class ResidualTest : public ::testing::TestWithParam<string>
{
protected:
    const string name = GetParam();
    const int n_points = 500;
    const real_t residual_tol = 1e-11;
    const real_t small_value_tol = 1e-300;

    void SetUp() override {}
    void TearDown() override {}

    void assert_small_residuals_(const string &intro_string, const vector<real_t> &Ra, const vector<real_t> &Rr)
    {
        const auto max_iter_abs = std::max_element(Ra.begin(), Ra.end());
        const auto max_iter_rel = std::max_element(Rr.begin(), Rr.end());
        const int max_index_abs = max_iter_abs - Ra.begin();
        const int max_index_rel = max_iter_rel - Rr.begin();
        const real_t max_value_abs = *max_iter_abs;
        const real_t max_value_rel = *max_iter_rel;

        // print and assert result
        std::cout << std::scientific << std::setprecision(2);
        std::cout << intro_string << ": max residuals of (";
        std::cout << max_index_abs << ", " << max_value_abs << ") and (";
        std::cout << max_index_rel << ", " << max_value_rel << ")" << std::endl;
        EXPECT_LT(max_value_abs, residual_tol);
    }
};

//-----------------------------------------------------------------
// Tests for Split Stencils (order 2, 4, 6, 8)
//-----------------------------------------------------------------

using LGFOneUnbounded::lgf2_small_cutoff_;
using LGFOneUnbounded::lgf4_small_cutoff_;
using LGFOneUnbounded::lgf6_small_cutoff_;
using LGFOneUnbounded::lgf8_small_cutoff_;

using LGFOneUnbounded::lgf4_repeat_cutoff_;
using LGFOneUnbounded::lgf8_repeat_cutoff_;

const real_t c4 = LGFOneUnbounded::lgf4_c_star_;
const real_t c8 = LGFOneUnbounded::lgf8_c_star_;

/**
 * @brief test values for each stencil
 */
const std::map<string, vector<real_t>> split_stencil_values = {
    {"LGF2", {0, 1e-8, 1e-6, 1e-4, 0.99 * lgf2_small_cutoff_, 0.01 * lgf2_small_cutoff_, 1e-2, 1e-1, 1.0, 2.0, 4.0, 6.0, 8.0}},
    {"LGF4",
     {0,
      1e-8, 1e-6, 1e-4, 0.99 * lgf4_small_cutoff_, 1.01 * lgf4_small_cutoff_,
      1e-2, 1.0, 2.0, 2.9,
      c4 - 1.01 * lgf4_repeat_cutoff_, c4 - 0.99 * lgf4_repeat_cutoff_, c4 - 1e-6, c4, c4 + 1e-6, c4 + 0.99 * lgf4_repeat_cutoff_, c4 + 1.01 * lgf4_repeat_cutoff_,
      4.0, 6.0, 8.0, 10.0, 12.0}},
    {"LGF6", {0, 1e-8, 1e-6, 1e-4, 0.99 * lgf6_small_cutoff_, 1.01 * lgf6_small_cutoff_, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0}},
    {"LGF8", {0, 1e-8, 1e-6, 1e-4, 0.99 * lgf8_small_cutoff_, 1.01 * lgf8_small_cutoff_, 1.0, 2.0, 3.0, c8 - 1.01 * lgf8_repeat_cutoff_, c8 - 0.99 * lgf8_repeat_cutoff_, c8 - 1e-6, c8, c8 + 1e-6, c8 + 0.99 * lgf8_repeat_cutoff_, c8 + 1.01 * lgf8_repeat_cutoff_, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0}}};

class SplitResidualTest : public ResidualTest {};

TEST_P(SplitResidualTest, split_residuals)
{
    const SplitStencilInfo stencil = split_stencil_lookup.at(name);
    const vector<real_t> c_values = split_stencil_values.at(name);

    for (const auto &c : c_values)
    {

        vector<real_t> G = fill_split_lgf_values(stencil, c, n_points); //<! LGF values for -w, ..., n+w-1
        vector<real_t> Ra(n_points);                           //<! absolute residual for 0, ..., n-1
        vector<real_t> Rr(n_points);                           //<! relative residual for 0, ..., n-1
        calculate_split_residuals(stencil, c, n_points, G, Ra, Rr, small_value_tol);

        const string test_name = name + ", c = " + to_string(c);
        assert_small_residuals_(test_name, Ra, Rr);
    }
}

INSTANTIATE_TEST_SUITE_P(SplitResidualTest, SplitResidualTest, testing::Values("LGF2", "LGF4", "LGF6", "LGF8"));

//-----------------------------------------------------------------
// Tests for Mehrstellen Stencils (order 4, 6, both LHS and Full)
//-----------------------------------------------------------------

vector<array<real_t, 2>> mehrstellen_k_values = {
    {0, 0}, 
    {0.0005, -0.0005}, 
    {0.005, 0.01},
    {0.01, -0.03},
    {0.2, 0.2}, 
    {0.5, 0.4}, 
    {-1.5, 2.0}, 
    {2. * M_PI / 3., 2. * M_PI / 3.}, 
    {M_PI, M_PI / 3.},
    {M_PI, M_PI / 2.},
    {M_PI, M_PI}
};

class MehrstellenResidualTest : public ResidualTest {};

TEST_P(MehrstellenResidualTest, mehrstellen_residuals)
{
    const MehrstellenStencilInfo stencil = mehrstellen_stencil_lookup.at(name);

    for (const auto &k : mehrstellen_k_values)
    {
        vector<real_t> G = fill_mehrstellen_lgf_values(stencil, k, n_points); //<! LGF values for -w, ..., n+w-1
        vector<real_t> Ra(n_points);                                 //<! absolute residual for 0, ..., n-1
        vector<real_t> Rr(n_points);                                 //<! relative residual for 0, ..., n-1
        
        calculate_mehrstellen_residuals(stencil, k, n_points, G, Ra, Rr, small_value_tol);
        const string test_name = name + ", k = [" + to_string(k[0]) + ", " + to_string(k[1]) + "]";
        assert_small_residuals_(test_name, Ra, Rr);
    }
}

INSTANTIATE_TEST_SUITE_P(MehrstellenResidualTest, MehrstellenResidualTest, testing::Values("MEH4_LEFT", "MEH6_LEFT", "MEH4_FULL", "MEH6_FULL"));