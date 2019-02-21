#include <iostream>
#include <iomanip>
#include <vector>
#include <array>
#include <math.h>
#include <algorithm>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "mkl.h"
extern "C" {
  #include "spglib.h"
  #include "fftw3.h"
}

namespace constants {
	const double TWOPI = 2*std::acos(-1);
	const double ZERO = 1e-10;
	const double TOL = 1e-4;
	const double SYMPREC = 1e-4;
}

#ifdef RUN_SKW_TEST
#define SKW_TEST 1
#else
#define SKW_TEST 0
#endif

#ifdef DEBUG_FFT
#define DEBUG_DUMP_FFT(x) x
#else
#define DEBUG_DUMP_FFT(x)
#endif

#ifdef DEBUG_RHO
#define DEBUG_DUMP_RHO(x) x
#else
#define DEBUG_DUMP_RHO(x)
#endif

#ifdef DEBUG_H
#define DEBUG_DUMP_H(x) x
#else
#define DEBUG_DUMP_H(x)
#endif

#ifdef DEBUG_SYMOPS
#define DEBUG_DUMP_SYMOPS(x) x
#else
#define DEBUG_DUMP_SYMOPS(x)
#endif

#ifdef DEBUG_LATPOINTS
#define DEBUG_DUMP_LATPOINTS(x) x
#else
#define DEBUG_DUMP_LATPOINTS(x)
#endif

#ifdef DEBUG_S
#define DEBUG_DUMP_S(x) x
#else
#define DEBUG_DUMP_S(x)
#endif

#ifdef INFO
#define INFO_DUMP(x) x
#else
#define INFO_DUMP(x)
#endif

void fetch_sym_ops(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<int> &, std::vector<std::vector<std::vector<int> > > &);

void specify_sym_ops(std::vector<std::vector<std::vector<int> > > &);

void generate_lat_vec(std::vector<std::vector<double> > &, std::vector<std::vector<std::vector<int> > > &, double &, std::vector<int> &, std::vector<std::vector<int> > &, std::vector<std::vector<std::vector<int> > > &, std::vector<double> &);

void locate_symmetry_vectors(std::vector<int> &, std::vector<std::vector<std::vector<int> > > &, std::vector<std::vector<int> > &);

void generate_rho(std::vector<std::vector<int> > &, std::vector<double> &, std::vector<double> &, std::vector<double> &);

void set_c(std::vector<double> &);

void generate_s(std::vector<std::vector<int> > &, std::vector<std::vector<double> > &, std::vector<std::vector<std::vector<int> > > &, std::vector<std::vector<double> > &);

void generate_h_and_ssr(std::vector<std::vector<double> > &, std::vector<double> &, std::vector<double> &, std::vector<std::vector<double> > &ssr);

void generate_energies_diff(std::vector<std::vector<double> > &, std::vector<double> &);

int generate_lambdas(std::vector<std::vector<double> > &, std::vector<double> &, std::vector<double> &);

void generate_epsilons(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<double> &, std::vector<std::vector<double> > &);

void check_interpolated_energies(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

void interpolate(std::vector<std::vector<double> > &, std::vector<std::vector<std::vector<int> > > &, std::vector<std::vector<double> > &, std::vector<std::vector<double> > &, std::vector<std::vector<std::vector<double> > > &, std::vector<std::vector<double> > &, std::vector<int> &);

void test_skw();

void transpose_matrix_3x3(std::vector<std::vector<double> > &, std::vector<std::vector<double> > &);

void multiply_vector_matrix_3x3_d(std::vector<std::vector<double> > &, std::vector<double> &,std::vector<double> & );

void multiply_vector_matrix_3x3_i(std::vector<std::vector<int> > &, std::vector<int> &,std::vector<int> & );

void multiply_matrix_matrix_3x3_i(std::vector<std::vector<int> > &, std::vector<std::vector<in\
t> > &, std::vector<std::vector<int> > &);

double dot_vectors(std::vector<double> &, std::vector<int> &);

double length_of_vector_d(std::vector<double> &);

double length_of_vector_i(std::vector<int> &);

void flatten_nxn_matrix(std::vector<std::vector<double> > &, std::vector<double> &, int);

void flatten_nxn_matrix_pointer(std::vector<std::vector<double> > &, double *, int);

template <typename T> std::vector<size_t> sort_indexes(const std::vector<T> &v) {
  std::vector<size_t> index(v.size());
  for (size_t i=0;i!=index.size();++i) {
    index[i]=i;
  }
  std::sort(index.begin(),index.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return index;
}
