#include "lgf_one_unbounded.hpp"
#include <string>
#include <vector>
#include <map>

using LGFOneUnbounded::real_t;

/**
 * @brief information for a general dimension-split stencil
 */
struct SplitStencilInfo
{
    const int width = -1;             //!< width of the stencil
    real_t (*lgf)(int, real_t);       //!< function pointer to lgf(n, c)
    real_t (*sigma)(real_t);          //!< function pointer to sigma(k)
    const std::vector<real_t> coeffs = {}; //!< coefficients of the stencil
};

/**
 * @brief information for split stencils of order 2, 4, 6, and 8
 */
const std::map<std::string, SplitStencilInfo> split_stencil_lookup = {
    {"LGF2", {1, &LGFOneUnbounded::lgf2, &LGFOneUnbounded::lgf2_symbol, {2., -1.}}},
    {"LGF4", {2, &LGFOneUnbounded::lgf4, &LGFOneUnbounded::lgf4_symbol, {5. / 2., -4. / 3., 1. / 12.}}},
    {"LGF6", {3, &LGFOneUnbounded::lgf6, &LGFOneUnbounded::lgf6_symbol, {49. / 18., -3. / 2., 3. / 20., -1. / 90.}}},
    {"LGF8", {4, &LGFOneUnbounded::lgf8, &LGFOneUnbounded::lgf8_symbol, {205. / 72., -8. / 5., 1. / 5., -8. / 315., 1. / 560}}}};

/**
 * @brief information for a Mehrstellen stencil
 */
struct MehrstellenStencilInfo
{
    const int left_width = -1;                      //!< width of the left stencil
    const int right_width = -1;                     //!< width of the right stencil
    real_t (*lgf)(int, real_t, real_t);             //!< function pointer to lgf(n, k1, k2)
    void (*left_coeffs)(real_t *, const real_t[]);  //!< function to fill coefficients of the left stencil
    void (*right_coeffs)(real_t *, const real_t[]); //!< function to fill coefficients of the right stencil
};

void dummy_right_coeffs(real_t coeffs[1], const real_t y[2])
{
    coeffs[0] = 1.0;
}

/**
 * @brief information for Mehrstellen stencils of order 4 and 6
 */
const std::map<std::string, MehrstellenStencilInfo> mehrstellen_stencil_lookup = {
    {"MEH4_LEFT", {1, 0, &LGFOneUnbounded::meh4_left, &LGFOneUnbounded::meh4_left_coeffs, dummy_right_coeffs}},
    {"MEH4_FULL", {1, 1, &LGFOneUnbounded::meh4_full, &LGFOneUnbounded::meh4_left_coeffs, &LGFOneUnbounded::meh4_right_coeffs}},
    {"MEH6_LEFT", {1, 0, &LGFOneUnbounded::meh6_left, &LGFOneUnbounded::meh6_left_coeffs, dummy_right_coeffs}},
    {"MEH6_FULL", {1, 2, &LGFOneUnbounded::meh6_full, &LGFOneUnbounded::meh6_left_coeffs, &LGFOneUnbounded::meh6_right_coeffs}}};


std::vector<real_t> fill_split_lgf_values(const SplitStencilInfo &stencil, real_t c, int N)
{
    const int w = stencil.width;
    std::vector<real_t> G(N + 2*w);
    for (int i = 0; i < N + 2*w; i++)
    {
        G[i] = (*stencil.lgf)(i - w, c);
    }
    return G;
}

void calculate_split_residuals(const SplitStencilInfo &stencil, real_t c, int N, const std::vector<real_t> &G, std::vector<real_t> &Ra, std::vector<real_t> &Rr, real_t small_value_tol = 1e-300)
{
    const int w = stencil.width;
    for (int n = 0; n < N; n++)
    {
        const int i = n + w;
        double value = (stencil.coeffs[0] + c) * G[i];
        for (int j = 1; j <= w; j++)
        {
            value += stencil.coeffs[j] * (G[i - j] + G[i + j]);
        }
        const double exact = (n == 0) ? 1.0 : 0.0;
        Ra[n] = fabs(value - exact);
        Rr[n] = (fabs(G[i]) < small_value_tol) ? -1 : fabs(value - exact) / fabs(G[i]);
    }
}

std::vector<real_t> fill_mehrstellen_lgf_values(const MehrstellenStencilInfo &stencil, std::array<real_t, 2> k, int N)
{
    const int wl = stencil.left_width;
    std::vector<real_t> G(N + 2 * wl);
    for (int i = 0; i < N + 2*wl; i++)
    {
        const int n = i - wl;
        G[i] = (*stencil.lgf)(n, k[0], k[1]);
    }
    return G;
}

void calculate_mehrstellen_residuals(const MehrstellenStencilInfo &stencil, std::array<real_t, 2> k, int N, const std::vector<real_t> &G, std::vector<real_t> &Ra, std::vector<real_t> &Rr, real_t small_value_tol = 1e-300)
{
    const real_t y[2] = {pow(sin(k[0]/2), 2), pow(sin(k[1]/2), 2)};
    real_t left_coeffs[stencil.left_width + 1];
    stencil.left_coeffs(left_coeffs, y);
    real_t right_coeffs[stencil.right_width + 1];
    stencil.right_coeffs(right_coeffs, y);

    for (int n = 0; n < N; n++)
    {
        const int i = n + stencil.left_width;
        double value = left_coeffs[0] * G[i];
        for (int j = 1; j <= stencil.left_width; j++)
        {
            value += left_coeffs[j] * (G[i - j] + G[i + j]);
        }
        const double exact = (n <= stencil.right_width) ? right_coeffs[n] : 0.0;
        Ra[n] = fabs(value - exact);
        Rr[n] = (fabs(G[i]) < small_value_tol) ? -1 : fabs(value - exact) / fabs(G[i]);
    }
}