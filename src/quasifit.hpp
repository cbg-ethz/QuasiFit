#include <config.h>

/* C++ I/O */
#include <iostream>
#include <fstream>


/* C++ STL types */
#include <string>
#include <vector>
#include <map>
#include <iterator>


/* Standard types and limits */
#include <cstdio>
#include <cstdlib>
#include <limits>
#ifdef HAVE_CSTDINT	// Mac OS X does not have a separate cstdint header
	#include <cstdint>
#endif


/* Progress Bar */
#include "ezETAProgressBar.hpp"


/* Console handling */
#include <console_buf.hpp>
extern console_buf console_buffer;
extern std::ostream console;


/* GNU Scientific Library */
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>


/* Boost */
#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

// In case -fno-exceptions is desired
// only works with the very latest Boost releases
/*
   #define BOOST_NO_EXCEPTIONS
   namespace boost
   {
        void throw_exception(std::exception const&)
        {
                std::cout << "Unexpected exception!\n";
                exit(EXIT_FAILURE);
        }
   }
 */

/* Floating-point definitions */
typedef double NORMAL_DOUBLE;

// Set up the extended floating point type:
#if defined(HAVE_LONG_DOUBLE_PRECISION)
	// ordinary extended precision (long double; ~18 digits on x87)
	#define PRECISION_TYPE "Extended Double"
	typedef long double EXT_DOUBLE;

#elif defined(HAVE_QUAD_PRECISION)
	// Quad precision (~34 digits)
	#define PRECISION_TYPE "Quad Precision"
	#include <boost/multiprecision/float128.hpp>
	typedef boost::multiprecision::float128 EXT_DOUBLE;

#elif defined(HAVE_TRUEMULTI_PRECISION)
	// GMP arbitrary precision (~100 digits)
	#define PRECISION_TYPE "Arbitrary Precision"
	#include <boost/multiprecision/gmp.hpp>
	typedef boost::multiprecision::mpf_float_100 EXT_DOUBLE;

#else
	// ERROR
	#error You have not selected a floating-point type to use. This is probably because you have not run the configure script yet.
#endif


/* Eigen */
#include <Eigen/Dense>

typedef Eigen::Matrix<EXT_DOUBLE, Eigen::Dynamic, Eigen::Dynamic> MatrixED;
typedef Eigen::Matrix<EXT_DOUBLE, Eigen::Dynamic, 1> VectorED;

typedef Eigen::Matrix<NORMAL_DOUBLE, Eigen::Dynamic, Eigen::Dynamic> MatrixND;
typedef Eigen::Matrix<NORMAL_DOUBLE, Eigen::Dynamic, 1> VectorND;

typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> MatrixID;
typedef Eigen::Matrix<uint32_t, Eigen::Dynamic, 1> VectorID;

// Solver to use:
typedef Eigen::PartialPivLU<MatrixED> Solver;	// Satisfactory, should work in principle
// typedef Eigen::FullPivLU<MatrixED> Solver; // Slowest, but always applicable


// structs:
struct stats
{
	uint32_t count;
	NORMAL_DOUBLE freq;

	stats() : count(1) {}
};

typedef std::vector<std::string> seq_cont;
typedef std::vector<uint32_t> seq_inds;

// variables:
#include <vars.hpp>

/* functions */
// Misc functions:
uint64_t random_seed();
uint32_t hamming_distance(const std::string& A, const std::string& B);

void about(const char* program_name);
void counter_display();
void parse_arguments(int argc, char** argv);

// I/O functions:
void load_inputFile(seq_cont& sequences);
void write_output(const seq_cont& sequences);
void r_loader_from_file();

// MCMC functions:
void population_initialiser();

void generate_sample(gsl_rng* r, VectorED& new_sample, const uint32_t chain_id);
void convert_from_R_to_P(const VectorED& R_vector, VectorED& P_vector_unnorm, VectorED& P_vector);
void convert_from_P_to_R(const VectorED& P_vector, VectorED& R_vector);

void probability_sampler(boost::barrier& syn_barrier, uint32_t thread_no);
extern void (* MCMC)(boost::barrier& syn_barrier, uint32_t thread_no);

EXT_DOUBLE convert_from_probability_to_fitness__manifold(const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult);
extern EXT_DOUBLE (* fitness_space) (const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult);

/*
   inline void convert_from_M_to_S(const VectorED& M_vector, VectorED& S_vector)
   {
        S_vector = M_vector - VectorED::Ones(DIM);
   }
 */