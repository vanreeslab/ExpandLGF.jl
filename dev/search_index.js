var documenterSearchIndex = {"docs":
[{"location":"examples.html#Generate-C-code-for-far-field-LGF-Expansion","page":"Examples","title":"Generate C code for far-field LGF Expansion","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using ExpandLGF\nusing Symbolics\n\nstencil = Mehrstellen6()\nterms = 3\nexpansion_coefficients = ExpandLGF.far_field_expansion_coeffs(stencil, terms)\ncode = ExpandLGF.far_field_3d_expansion(expansion_coefficients, Symbolics.CTarget())","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The block above generates C code that evaluates a far-field expansion for the sixth-order Mehrstellen stencil. The length of the expansion is determined by terms, which indicates the number of terms m in the summation","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"G(n) = frac14pi bmn + frac1pisum_k=1^m fracg_k(bmn)bmn^6k+1 + mathcalOleft( bmn^-(2m + 3) right)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The expansion computed above will contain the analytical Green's function and the corrections g_1(bmn), g_2(bmn), and g_3(bmn), which provides an accuracy of mathcalO(n^-9). The output code is the string","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"#include <math.h>\nvoid lgf(double* G, const double n1, const double n2, const double n3, const double n) {\n  G[0] = (-0.018229166666666668 * pow(n1, 12) + -0.018229166666666668 * pow(n2, 12) + -0.018229166666666668 * pow(n3, 12) + -0.5625 * (n1 * n1) * pow(n2, 10) + 8.088541666666666 * pow(n1, 6) * pow(n2, 6) + 3.5 * pow(n1, 4) * pow(n2, 8) + 8.088541666666666 * pow(n1, 6) * pow(n3, 6) + -0.5625 * (n1 * n1) * pow(n3, 10) + 3.5 * pow(n1, 8) * pow(n2, 4) + 3.5 * pow(n1, 8) * pow(n3, 4) + -0.5625 * pow(n1, 10) * (n2 * n2) + 3.5 * pow(n1, 4) * pow(n3, 8) + -0.5625 * pow(n1, 10) * (n3 * n3) + 3.5 * pow(n2, 4) * pow(n3, 8) + 3.5 * pow(n2, 8) * pow(n3, 4) + 8.088541666666666 * pow(n2, 6) * pow(n3, 6) + -0.5625 * pow(n2, 10) * (n3 * n3) + -0.5625 * (n2 * n2) * pow(n3, 10) + -15.875 * pow(n1, 6) * pow(n2, 4) * (n3 * n3) + -9.953125 * pow(n1, 8) * (n2 * n2) * (n3 * n3) + -15.875 * pow(n1, 4) * (n2 * n2) * pow(n3, 6) + -15.875 * pow(n1, 6) * (n2 * n2) * pow(n3, 4) + -33.1640625 * pow(n1, 4) * pow(n2, 4) * pow(n3, 4) + -15.875 * pow(n1, 4) * pow(n2, 6) * (n3 * n3) + -9.953125 * (n1 * n1) * (n2 * n2) * pow(n3, 8) + -15.875 * (n1 * n1) * pow(n2, 6) * pow(n3, 4) + -15.875 * (n1 * n1) * pow(n2, 4) * pow(n3, 6) + -9.953125 * (n1 * n1) * pow(n2, 8) * (n3 * n3)) / (M_PI * pow(n, 19)) + 1 / (12.566370614359172 * n * 1);\n}","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"Note that because the stencil is sixth order, g_1(bmn) = g_2(bmn) = 0. For lower order stencils theses terms are nonzero and will appear in the generated C code.","category":"page"},{"location":"examples.html#Precompute-a-block-of-2D-LGF-values","page":"Examples","title":"Precompute a block of 2D LGF values","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using ExpandLGF\nusing OffsetArrays\n\nN = 8            # Evaluate entries on [0, N - 1]^2\nnear_terms = 10  # number of terms in near-field expansions\nfar_terms = 8    # number of terms used in far-field expansion\nrtol = 0         # no relative tolerance on evaluations\natol = 1e-15     # absolute tolerance on evaluations\n\nstencil = StandardDifference2D(4)\nELGF = EvaluateLGF(stencil, near_terms, far_terms, rtol, atol)\n    \nlgf = OffsetArray(zeros(Float64, (N, N)), 0:N-1, 0:N-1)\nfor index in ExpandLGF.sorted_indices_2D(N - 1)\n    lgf[index...] = ELGF(index)\nend\n\nExpandLGF.symmetrize_block!(lgf)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"Output:","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"8×8 OffsetArray(::Matrix{Float64}, 0:7, 0:7) with eltype Float64 with indices 0:7×0:7:\n  0.0       -0.207425  -0.318808  …  -0.495317  -0.519862\n -0.207425  -0.269189  -0.338713     -0.497503  -0.521472\n -0.318808  -0.338713  -0.37609      -0.503719  -0.526114\n -0.384607  -0.393345  -0.414417     -0.513102  -0.533294\n -0.430685  -0.435572  -0.448589     -0.524612  -0.542364\n -0.466274  -0.469411  -0.478132  …  -0.537315  -0.552685\n -0.495317  -0.497503  -0.503719     -0.550506  -0.563713\n -0.519862  -0.521472  -0.526114     -0.563713  -0.575037","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The code above evaluates a block of LGF entries for the 2D dimension-split stencil of order four. To take advantage of symmetry, the sorted_indices_2D function limits the loop to indices (i j) with i le j. Once the computation is complete, the rest of the block is filled using the symmetry condition G(i j) = G(j i) in the call to ExpandLGF.symmetrize_block!. ","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"For a 3D stencil, 2D can be replaced with 3D in the code and the array lgf can be allocated with","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"lgf = OffsetArray(zeros(Float64, (N, N, N)), 0:N-1, 0:N-1, 0:N-1)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"In three dimensions, the symmetry condition reduces the number of evaluations from N^3 to roughly N^3  6.","category":"page"},{"location":"index.html#ExpandLGF.jl","page":"Overview","title":"ExpandLGF.jl","text":"","category":"section"},{"location":"index.html","page":"Overview","title":"Overview","text":"ExpandLGF.jl provides fast and accurate evaluations of the Lattice Green's Functions (LGFs) for high-order finite difference discretizations of Poisson's equation.","category":"page"},{"location":"index.html#Background","page":"Overview","title":"Background","text":"","category":"section"},{"location":"index.html","page":"Overview","title":"Overview","text":"Lattice Green's Functions (LGFs) are fundamental solutions to discretized linear operators that are frequently used in the numerical solution to elliptic PDEs. As such, they are a useful tool for solving Poisson's equation on domains with one or more free-space boundary conditions, which appear frequently in simulations of incompressible fluid flows (see ViscousFlow.jl for one such example). ","category":"page"},{"location":"index.html","page":"Overview","title":"Overview","text":"Generally LGFs are defined through a Fourier integral with a singular and oscillatory integrand. For the standard second-order centered discretization of the Laplacian, analytical techniqes can be used to greatly simplify the 2D or 3D integral into a 1D form which is amenable to numerical quadrature. This package uses a combination of symbolic computation and numerical quadrature to extend these simplifications to higher-order discretizations. ","category":"page"},{"location":"index.html","page":"Overview","title":"Overview","text":"For interactions between distant sources, a far-field asymptotic expansion of the LGF can bypass numerical integration altogether. This package implements an existing algorithm from Martinsson and Rodin to calculate these far-field expansions. The results are available either as compiled Julia functions or C code.","category":"page"},{"location":"index.html","page":"Overview","title":"Overview","text":"Finally, the LGF concept can be extended to domains with a combination of periodic and free-space boundary conditions. For domains with one unbounded dimension and one or more periodic dimensions, a C++ header that evaluates the LGFs of several common discretizations is included in the ccode folder.","category":"page"},{"location":"docstrings.html#Docstrings","page":"Docstrings","title":"Docstrings","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"CurrentModule = ExpandLGF","category":"page"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"There are two ways to use this package for LGF evaluation. ","category":"page"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"The high level interface performs LGF evaluations by automatically switching from quadrature to a series expansion whenever it is permissible within a given error tolerance. The user specifies the tolerance and the number of series terms used. This functionality is exported by the package.\nThe low level interface provides direct access to expansion coefficients, function generation, and error estimates for each type of series. This functionality is not exported.","category":"page"},{"location":"docstrings.html#High-level-interface","page":"Docstrings","title":"High level interface","text":"","category":"section"},{"location":"docstrings.html#Stencil-Types","page":"Docstrings","title":"Stencil Types","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"Stencil\nSplitStencil\nFullStencil\nMehrstellenStencil\nMehrstellen4\nMehrstellen6\nLeftMehrstellen4\nLeftMehrstellen6\nStandardDifference2D\nStandardDifference3D","category":"page"},{"location":"docstrings.html#ExpandLGF.Stencil","page":"Docstrings","title":"ExpandLGF.Stencil","text":"An abstract type for all finite difference stencils\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.SplitStencil","page":"Docstrings","title":"ExpandLGF.SplitStencil","text":"A dimension-split finite difference stencil\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.FullStencil","page":"Docstrings","title":"ExpandLGF.FullStencil","text":"A general symmetric finite difference stencil in 2D or 3D\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.MehrstellenStencil","page":"Docstrings","title":"ExpandLGF.MehrstellenStencil","text":"A Mehrstellen stencil with both left and right finite difference stencils\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.Mehrstellen4","page":"Docstrings","title":"ExpandLGF.Mehrstellen4","text":"Return the 4th order 3D Mehrstellen stencil\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.Mehrstellen6","page":"Docstrings","title":"ExpandLGF.Mehrstellen6","text":"Return the 6th order 3D Mehrstellen stencil\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.LeftMehrstellen4","page":"Docstrings","title":"ExpandLGF.LeftMehrstellen4","text":"Return the stencil for the 4th order left-hand side Mehrstellen operator\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.LeftMehrstellen6","page":"Docstrings","title":"ExpandLGF.LeftMehrstellen6","text":"Return the stencil for the 6th order left-hand side Mehrstellen operator\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.StandardDifference2D","page":"Docstrings","title":"ExpandLGF.StandardDifference2D","text":"Return a 2D dimension-split stencil of given order and minimal width\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.StandardDifference3D","page":"Docstrings","title":"ExpandLGF.StandardDifference3D","text":"Return a 3D dimension-split stencil of given order and minimal width\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#LGF-Evaluation","page":"Docstrings","title":"LGF Evaluation","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"EvaluateLGF\nEvaluateLGF(n)\nNearFieldEvaluator\nNearFieldSplit\nNearFieldFull\nNearField\nFarFieldExpansion\nFarField","category":"page"},{"location":"docstrings.html#ExpandLGF.EvaluateLGF","page":"Docstrings","title":"ExpandLGF.EvaluateLGF","text":"A callable struct wrapping near-field evaluation, far-field evaluation, and a cutoff between the two. Call via (::EvaluateLGF)(n) for a 2D or 3D integer vector n. \n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.EvaluateLGF-Tuple{Any}","page":"Docstrings","title":"ExpandLGF.EvaluateLGF","text":"EvaluateLGF(stencil::SplitStencil, near_terms, far_terms, rtol, atol)\n\nReturn a callable struct that evaluates the LGF by switching between near-field and far-field evaluation strategies.\n\nArguments\n\nstencil : the split finite difference stencil\nnear_terms : number of inner and outer integral expansion terms.\nfar_terms : number of far-field expansion terms.\nrtol, atol : relative and absolute error tolerance for the evaluation.\n\n\n\n\n\nEvaluateLGF(stencil::Union{FullStencil, MehrstellenStencil}, far_terms, rtol, atol, maxevals)\n\nReturn a callable struct that evaluates the LGF for a FullStencil or MehrstellenStencil, switching between near-field or far-field evaluation strategies.\n\nArguments\n\nstencil : the finite difference stencil\nfar_terms : number of far-field expansion terms.\nrtol, atol : relative and absolute error tolerance for the evaluation.\nmaxevals : the max number of symbol evaluations used to compute near-field LGF values\n\n\n\n\n\n","category":"method"},{"location":"docstrings.html#ExpandLGF.NearFieldEvaluator","page":"Docstrings","title":"ExpandLGF.NearFieldEvaluator","text":"A callable structure that evaluates near-field values through quadrature. Call with (::NearFieldEvaluator)(n) for an integer vector n. \n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.NearFieldSplit","page":"Docstrings","title":"ExpandLGF.NearFieldSplit","text":"A callable structure that evaluates near-field values for a dimension-split stencil\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.NearFieldFull","page":"Docstrings","title":"ExpandLGF.NearFieldFull","text":"A callable structure that evaluates near-field values via quadrature for a full or Mehrstellen stencil\n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.NearField","page":"Docstrings","title":"ExpandLGF.NearField","text":"NearField(stencil::SplitStencil, nmax, terms, rtol, atol)\n\nReturn a callable struct which evaluates the LGF for a SplitStencil for n_infty le n_mathrmmax with prescribed error tolerance. Uses a combination of direct quadrature and series expansions.\n\nArguments\n\nstencil : the finite difference stencil\nnmax : the largest component of n that will be evaluated\nterms : number of terms used in integral expansions. Affects performance.\nrtol, atol : relative and absolute error tolerance for LGF values\n\n\n\n\n\nNearField(stencil::Stencil, rtol, atol, maxevals)\n\nReturn a callable struct which evaluates the LGF for any Stencil with prescribed error tolerance. Uses direct quadrature.\n\nArguments\n\nstencil : the finite difference stencil\nnmax : the largest component of n that will be evaluated\nrtol, atol : relative and absolute error tolerance for LGF values\nmaxevals : the max number of symbol evaluations used to compute each LGF value\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.FarFieldExpansion","page":"Docstrings","title":"ExpandLGF.FarFieldExpansion","text":"A callable struct for far field evaluation. Call via (::FarFieldExpansion)(n) for a 2D or 3D integer vector n. \n\n\n\n\n\n","category":"type"},{"location":"docstrings.html#ExpandLGF.FarField","page":"Docstrings","title":"ExpandLGF.FarField","text":"far_field FarField(stencil::Stencil, nterms)\n\nEvaluate the LGF for stencil via a far-field series expansion with ntemrs terms.\n\n\n\n\n\nfar_field, nmin = FarField(stencil::Stencil, nterms, rtol, atol)\n\nEvaluate the LGF for n_2 ge n_mathrmmin via series expansion with prescribed error tolerance. Returns a callable struct and the threshold nmin beyond which the expansion achieves the given tolerance.\n\nArguments\n\nstencil - the finite difference stencil\ndimensions : either 2 or 3 to indicate a 2D or 3D LGF\nnterms - number of terms in the expansion\nrtol, atol - relative and absolute error tolerance on the evaluation (pass 0 to ignore either).\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#Low-level-interface","page":"Docstrings","title":"Low level interface","text":"","category":"section"},{"location":"docstrings.html#Expansion-coefficients","page":"Docstrings","title":"Expansion coefficients","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"inner_integral_expansion_coeffs\nouter_integral_expansion_coeffs\nouter_integral_2d_expansion_coeffs\nouter_integral_3d_expansion_coeffs\ngraded_series_expansion\nfar_field_expansion_coeffs\nmartinsson_rodin_a2q","category":"page"},{"location":"docstrings.html#ExpandLGF.inner_integral_expansion_coeffs","page":"Docstrings","title":"ExpandLGF.inner_integral_expansion_coeffs","text":"inner_integral_expansion_coeffs(stencil::SplitStencil, n, m)\n\nReturn an m-term expansion of the inner integral for fixed n and large t. The results is an array b of univariate polynomials in n so that \n\nI_sigma n(t) = frac1sqrt4pi t sum_k=0^m b_k(n) t^-k + mathcalOleft( t^-frac2m+32 right)\n\nArguments\n\nstencil : the finite difference stencil\nm : the number of terms in the expansion\nPR : a univariate PolynomialRing\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_expansion_coeffs","page":"Docstrings","title":"ExpandLGF.outer_integral_expansion_coeffs","text":"outer_integral_expansion_coeffs(stencil::SplitStencil{3, T}, n, m) where T\n\nReturn an m-term expansion of the 3D outer integral for large T. The results is an array g of mulitvariate polynomials in n so that \n\nint_T^infty I_sigma n_1(t) I_sigma n_2(t) I_sigma n_3(t) mathrmdt approx frac1sqrt16pi^3T sum_k=0^m g_k(bmn) T^-k + mathcalOleft( T^-frac2m+32 right)\n\nArguments\n\nstencil : the finite difference stencil\nn : a vector of 3 generators from a multivariate PolynomialRing\nm : the number of terms in the expansion\n\n\n\n\n\nouter_integral_expansion_coeffs(stencil::SplitStencil{2, T}, n, m) where T\n\nReturn an m-term expansion of the 2D outer integral for large T. The results is an array g of mulitvariate polynomials in n so that \n\nint_T^infty I_sigma n_1(t) I_sigma n_2(t) - I_sigma 0(t)^2 mathrmdt approx frac14pi sum_k=1^m g_k(bmn) T^-k + mathcalOleft( T^-(m+1) right)\n\nArguments\n\nstencil : the finite difference stencil\nn : a vector of 2 generators from a multivariate PolynomialRing\nm : the number of terms in the expansion\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_2d_expansion_coeffs","page":"Docstrings","title":"ExpandLGF.outer_integral_2d_expansion_coeffs","text":"outer_integral_2d_expansion_coeffs(b, PR::MPolyRing)\n\nReturn a power series expansion of the 2D outer integral based on an expansion of the inner integral.\n\nArguments\n\nb : inner expansion, the result of inner_integral_expansion_coeffs()\nPR : a multivariate PolynomialRing in two variables\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_3d_expansion_coeffs","page":"Docstrings","title":"ExpandLGF.outer_integral_3d_expansion_coeffs","text":"outer_integral_3d_expansion_coeffs(b, n)\n\nReturn a power series expansion of the 3D outer integral based on an expansion of the inner integral.\n\nArguments\n\nb : inner expansion, the result of inner_integral_expansion_coeffs()\nPR : a multivariate PolynomialRing with 3 generators\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.graded_series_expansion","page":"Docstrings","title":"ExpandLGF.graded_series_expansion","text":"bj = graded_series_expansion(stencil::Stencil{T}, j_max::Int, PR::MPolyRing{T}) where T\n\nReturn a power series expansion of the Fourier symbol for a given stencil:\n\nsigma(bmk) = sum_j = 0^j_mathrmmax b_j(bmk) + mathcalOleft( bmk^2j_mathrmmax + 2 right)\n\nIn the output, b[j] contains a homogeneous multivariate polynomial of degree 2j.\n\nArguments\n\nstencil : the stencil defining the symbol\nj_max : the number of valid entries in the output bj\nPR : a multivariate PolynomialRing defining the type of the expansion\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.far_field_expansion_coeffs","page":"Docstrings","title":"ExpandLGF.far_field_expansion_coeffs","text":"far_field_expansion_coeffs(stencil::Stencil{3, T}, m, PR::MPolyRing{T}) where T\n\nCalculates an m-term far-field expansion of the LGF using an algorithm from Martinsson and Rodin. The result is an array g of multivariate polynomials in n so that\n\nG(n) = frac14pi bmn + frac1pisum_k=1^m fracg_k(bmn)bmn^6k+1 + mathcalOleft( bmn^-(2m + 3) right)\n\nNote that once m exceeds 5 or 6, each additional term will roughly double the runtime of this function.\n\nArguments\n\nstencil : the 3D finite difference stencil\nm : number of correction terms returned\nPR : a multivariate PolynomialRing that defines the type of the output\n\n\n\n\n\nfar_field_expansion_coeffs(stencil::Stencil{2, T}, m, PR::MPolyRing{T}) where T\n\nCalculates an m-term far-field expansion of the 2D LGF using an algorithm from Martinsson and Rodin. The result is an array g of multivariate polynomials in n so that\n\nG(bmn) = -frac12pileft(log bmn + gamma + fraclog 82right) + frac1pisum_k=1^m fracg_k(bmn)bmn^6j + mathcalOleft( bmn^-(m+1) right)\n\nNote that once m exceeds 5 or 6, each additional term will roughly double the runtime of this function.\n\nArguments\n\nstencil : the 2D finite difference stencil\nm : number of correction terms returned\nPR : a multivariate PolynomialRing that defines the type of the output\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.martinsson_rodin_a2q","page":"Docstrings","title":"ExpandLGF.martinsson_rodin_a2q","text":"a2q = martinsson_rodin_a2q(stencil::Stencil, q_max::Int, PR::MPolyRing)\n\nReturn the polynomials a_2q(bmk) with q <= q_max for a given stencil.\n\nSee Martinsson, Rodin. \"Asymptotic Expansions of Lattice Green's Functions\". JCP 2022\n\nArguments\n\ncffs : the finite difference coefficients\nq_max : the first q_max polynomials are returned\nPR : a vector of 3 generators from a multivariate PolynomialRing\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#Callable-expansions","page":"Docstrings","title":"Callable expansions","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"inner_integral_expansion\nouter_integral_expansion\nouter_integral_2d_expansion\nouter_integral_3d_expansion\nfar_field_expansion\nfar_field_2d_expansion\nfar_field_3d_expansion","category":"page"},{"location":"docstrings.html#ExpandLGF.inner_integral_expansion","page":"Docstrings","title":"ExpandLGF.inner_integral_expansion","text":"inner_integral_expansion(stencil::SplitStencil, terms)\n\nReturn a function which evaluates I_sigman(t) for large t via a series expansion.\n\nArguments\n\nstencil : coefficients of the FD scheme\nterms : the number of terms in the expansion\n\nReturn\n\nout_fun : a scalar function with signature val = out_fun(n, t)\n\n\n\n\n\ninner_integral_expansion(b)\n\nReturn a function which evaluates I_sigman(t) for large t via a series expansion.\n\nArguments\n\nb : expansion coefficients, returned by inner_integral_expansion_coeffs\n\nReturn\n\nout_fun : a scalar function with signature val = out_fun(n, t)\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_expansion","page":"Docstrings","title":"ExpandLGF.outer_integral_expansion","text":"outer_integral_expansion(stencil::SplitStencil{N, T}, terms) where {N, T}\n\nReturn a Julia function which evaluates the outer integral over tinfty for large t via a series expansion.\n\nArguments\n\nstencil : coefficients of the FD scheme\nterms : the number of terms in the expansion\n\nReturn\n\nout_fun : a scalar function with signature val = out_fun(n, t) for an N-Vector n and scalar t.\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_2d_expansion","page":"Docstrings","title":"ExpandLGF.outer_integral_2d_expansion","text":"outer_integral_2d_expansion(g)\n\nReturn a Julia function which evaluates the 2D outer integral over Tinfty for large T via a series expansion.\n\nArguments\n\ng : expansion coefficients, the result of outer_integral_2d_expansion_coeffs()\n\nReturn\n\nout_fun : a scalar function with signature val = out_fun([n1, n2], T)\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_3d_expansion","page":"Docstrings","title":"ExpandLGF.outer_integral_3d_expansion","text":"outer_integral_3d_expansion(g)\n\nReturn a Julia function which evaluates the 3D outer integral over Tinfty for large T via a series expansion.\n\nArguments\n\ng : expansion coefficients, the result of outer_integral_3d_expansion_coeffs()\n\nReturn\n\nout_fun : a scalar function with signature val = out_fun([n1, n2, n3], T)\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.far_field_expansion","page":"Docstrings","title":"ExpandLGF.far_field_expansion","text":"far_field_expansion(stencil::Stencil{N, T}, nterms, target = Symbolics.JuliaTarget()) where {N, T}\n\nReturn a function which evaluates a far-field expansion of the LGF.\n\nWhen target = Symbolics.JuliaTarget() (default), outputs a function G = out_fun(n) for scalar G, N-vector n. When target = Symbolics.CTarget(), outputs a string containing a C function that evaluates the expansion.\n\nArguments\n\nstencil : coefficients of the finite difference stencil\nnterms : number of terms in the expansion\ntarget : choose either Symbolics.JuliaTarget()` or Symbolics.CTarget()\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.far_field_2d_expansion","page":"Docstrings","title":"ExpandLGF.far_field_2d_expansion","text":"far_field_2d_expansion(g, target::Symbolics.JuliaTarget)\n\nReturn a Julia function which evaluates a far-field expansion of the 2D LGF with coefficients in g. The output has signature G = out_fun(n) for scalar G, 2-vector n\n\nArguments\n\ng : expansion coefficients, returned by far_field_2d_expansion_coeffs\n\n\n\n\n\nfar_field_2d_expansion(g, target::Symbolics.CTarget)\n\nReturn a string containing C code that evaluates a far-field expansion of the 2D LGF.\n\nArguments\n\ng : expansion coefficients, returned by far_field_2d_expansion_coeffs\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.far_field_3d_expansion","page":"Docstrings","title":"ExpandLGF.far_field_3d_expansion","text":"far_field_3d_expansion(g, target::JuliaTarget)\n\nReturn a Julia function which evaluates a far-field expansion of the 3D LGF. The output has signature G = out_fun(n) for scalar G, 3-vector n\n\nArguments\n\ng : expansion coefficients, returned by far_field_3d_expansion_coeffs\n\n\n\n\n\nfar_field_3d_expansion(g, target::CTarget)\n\nReturn a string containing C code that evaluates a far-field expansion of the 3D LGF.\n\nArguments\n\ng : far-field expansion coefficients, returned by far_field_3d_expansion_coeffs\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#Error-estimates","page":"Docstrings","title":"Error estimates","text":"","category":"section"},{"location":"docstrings.html","page":"Docstrings","title":"Docstrings","text":"inner_integral_errors\nouter_integral_errors\nfar_field_errors","category":"page"},{"location":"docstrings.html#ExpandLGF.inner_integral_errors","page":"Docstrings","title":"ExpandLGF.inner_integral_errors","text":"tmax = inner_integral_errors(b, nmax; rtol=1e-16, atol=1e-16)\n\nCalculate the region of validity for an inner integral expansion.\n\nSpecifically, when n le n_mathrmmax and t ge t_mathrmmax, the last term in the  expansion is smaller than both atol and rtol * (first term).\n\nArguments\n\nb : vector of expansion coefficients (polynomials in n)\nnmax : the maximal value of n that will be used in the expanions\nrtol : relative tolerance on the last term in the series (pass 0 to ignore)\natol : absolute tolerance on the last term in the series (pass 0 to ignore)\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.outer_integral_errors","page":"Docstrings","title":"ExpandLGF.outer_integral_errors","text":"Tmax = outer_integral_errors(g, nmax; rtol=1e-16, atol=1e-16)\n\nCalculate the region of validity for an outer integral expansion. The dimension (2D or 3D) is inferred from the type of g.\n\nSpecifically, when n_infty  n_mathrmmax and T ge T_mathrmmax, the last term in the  expansion is smaller than both atol and rtol * (first term).\n\nArguments\n\ng : vector of expansion coefficients (polynomials in n1, n2, n3)\nnmax : the maximal value of the 3-vector n that will be used in the expanions\nrtol : relative tolerance on the last term in the series (pass 0 to ignore)\natol : absolute tolerance on the last term in the series (pass 0 to ignore)\n\n\n\n\n\n","category":"function"},{"location":"docstrings.html#ExpandLGF.far_field_errors","page":"Docstrings","title":"ExpandLGF.far_field_errors","text":"nmin = far_field_errors(g; rtol=1e-16, atol=1e-16, search=100)\n\nCalculate the region of validity for a far-field expansion. The dimension (2D or 3D) is inferred from the type of g.\n\nSpecifically, when n_2 ge n_mathrmmin, the last term in the expansion is smaller than both  atol and rtol * (first term).\n\nArguments\n\ng : vector of far-field expansion coefficients (polynomials in n1, n2, n3)\nrtol : relative tolerance on the last term in the series (pass 0 to ignore)\natol : absolute tolerance on the last term in the series (pass 0 to ignore)\nsearch : finding nmin requires a bisection search over n. The search arg provides an upper limit.\n\n\n\n\n\n","category":"function"}]
}