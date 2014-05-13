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
#ifdef HAVE_CSTDINT // Mac OS X does not have a separate cstdint header
	#include <cstdint>
#endif


/* Progress Bar */
#include "ezETAProgressBar.hpp"


/* Console handling */
#include <console_buf.hpp>
extern console_buf console_buffer;
extern std::ostream console;

/* Command-line parsing */
#include <getopt.h>
static struct option long_options[] =
{
	{0, 0, 0, 0}
};


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
	namespace boost {
		void throw_exception(std::exception const&) {
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

typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic> MatrixID;
typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> VectorID;

// Solver to use:
typedef Eigen::PartialPivLU<MatrixED> Solver; // Satisfactory, should work in principle
//typedef Eigen::FullPivLU<MatrixED> Solver; // Slowest, but always applicable


/* Constants */
extern const uint64_t alph_card;
extern const EXT_DOUBLE global_m; // global HIV mutation rate
extern const EXT_DOUBLE base_u; // A -> T mutation rate

extern const char delim;


/* Command-line options */
extern uint64_t N;
extern uint64_t T;
extern uint64_t C;
extern uint64_t s;

extern uint64_t verbosity_level;

extern bool use_fitness_sampler;
extern bool use_simplex_space;

extern bool no_headers;
extern bool skip_writing;
extern bool omit_empties;
extern bool help_flag;
extern bool load_initial_r_from_file;
extern bool randomise_leave_out;
extern bool MLE_exists;

extern NORMAL_DOUBLE GAMMA;
extern NORMAL_DOUBLE B;

extern std::string inputFile;


/* Global variables */
extern uint64_t DIM;
extern uint64_t Nreads;
extern uint64_t L;
extern VectorID counters;
extern VectorID acc;
extern VectorID shuffled_index;

extern uint32_t leave_id;

extern boost::mutex mutex;

extern VectorID Data;

extern Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;
extern MatrixED Q_Random;

extern MatrixED QT;
extern MatrixED sQTinv;
extern MatrixED fQTinv;

extern MatrixED matrix_f;
extern MatrixED matrix_m;
extern MatrixED matrix_s;
extern MatrixED matrix_p;
extern MatrixED matrix_r;
extern MatrixED matrix_ln;

extern MatrixED population_f;
extern MatrixED population_m;
extern MatrixED population_s;
extern MatrixED population_p;
extern MatrixED population_r;
extern MatrixED population_r_new;
extern MatrixED population_ln;
extern MatrixED population_ln_new;


// MCMC options
extern uint64_t NChains_per_thread;
extern uint64_t no_samples;


// I/O options, defined in io.cpp
extern std::ostream console;

extern std::string inputFile_initial;
extern std::string output_f_File;
extern std::string output_m_File;
extern std::string output_s_File;
extern std::string output_p_File;
extern std::string output_r_File;
extern std::string output_ln_File;


// structs:
struct stats
{
	uint64_t count;
	NORMAL_DOUBLE freq;
	bool write_out;
	
	stats() : count(1), freq(0), write_out(true) {}
};

typedef std::vector<std::string> seq_cont;
typedef std::vector<uint32_t> seq_inds;

// functions:
uint64_t hamming_distance(const std::string& A, const std::string& B);

void about(const char* program_name);

void load_inputFile(seq_cont& sequences, seq_inds& indices, VectorED& p_MLE);

uint64_t random_seed();

void counter_display();

void generate_sample(gsl_rng* r, VectorED& new_sample, const uint64_t chain_id);

EXT_DOUBLE convert_from_probability_to_fitness__manifold(const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, VectorED& M_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult);

inline void convert_from_M_to_S(const VectorED& M_vector, VectorED& S_vector)
{
	S_vector = M_vector - VectorED::Ones(DIM);
}

extern EXT_DOUBLE (*fitness_space) (const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, VectorED& M_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult);

void convert_from_R_to_P(const VectorED& R_vector, VectorED& P_vector_unnorm, VectorED& P_vector);

void convert_from_P_to_R(const VectorED& P_vector, VectorED& R_vector);

void probability_sampler(boost::barrier& syn_barrier, uint64_t thread_no);

extern void (*MCMC)(boost::barrier& syn_barrier, uint64_t thread_no);// = probability_sampler;

void r_loader_from_file();

void population_initialiser(uint64_t thread_no);

void parse_arguments(int argc, char** argv);
