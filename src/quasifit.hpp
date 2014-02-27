#include <iostream>
#include <iomanip>
#include <fstream>

#include <string>
#include <vector>
#include <map>

#include <iterator>

#include <cstdio>
#include <cstdlib>

#ifdef __linux__
#include <cstdint>
#endif

//#include <cmath>
#include <limits>

#include <getopt.h>
static struct option long_options[] =
{
    {"sF",      no_argument, 0, 1000},
    {"simplex", no_argument, 0, 1001},
    {"nh",      no_argument, 0, 1010},
    {"nR",      no_argument, 0, 1020},
    {"initial", no_argument, 0, 1030},
    {0, 0, 0, 0}
};

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>

#include <boost/chrono.hpp>
#include <boost/chrono/process_cpu_clocks.hpp>
#include <boost/algorithm/string.hpp>

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
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>

// typedefs:
typedef double NORMAL_DOUBLE;

// Set up the extended floating point type:
#if defined(HAVE_QUAD_PRECISION)
    // Quad precision
#define PRECISION_TYPE Quad Precision

typedef _Quad EXT_DOUBLE;

namespace std {
    inline _Quad abs (const _Quad& q) { return __fabsq(q); }
    inline _Quad sqrt(const _Quad& q) { return __sqrtq(q); }
    inline _Quad log (const _Quad& q) { return __logq(q) ; }
    inline _Quad ceil(const _Quad& q) { return __ceilq(q); }

    ostream& operator<<( ostream& output, const EXT_DOUBLE& q )
    {
        output << static_cast<double>(q);
        return output;
    }
}

#elif defined(HAVE_MULTI_PRECISION)
    // MPFR arbitrary precision
#define PRECISION_TYPE Arbitrary Precision
    static_assert(0, "MULTI!");
#else
    // ordinary long double
#define PRECISION_TYPE Extended Double
    typedef long double EXT_DOUBLE;
#endif

#include <Eigen/Dense>

typedef Eigen::Matrix<EXT_DOUBLE, Eigen::Dynamic, Eigen::Dynamic> MatrixED;
typedef Eigen::Matrix<EXT_DOUBLE, Eigen::Dynamic, 1> VectorED;

typedef Eigen::Matrix<NORMAL_DOUBLE, Eigen::Dynamic, Eigen::Dynamic> MatrixND;
typedef Eigen::Matrix<NORMAL_DOUBLE, Eigen::Dynamic, 1> VectorND;

typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, Eigen::Dynamic> MatrixID;
typedef Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> VectorID;

/* Solver to use: */
typedef Eigen::PartialPivLU<MatrixED> Solver; // Satisfactory, should always work
//typedef Eigen::FullPivLU<MatrixED> Solver; // Slowest, but always applicable
//typedef Eigen::LLT<MatrixED> CholeskySolver; // Fastest, but only works on full sequence space

// Constants:
const uint64_t alph_card = 4;
const EXT_DOUBLE global_m = 3E-5; // global HIV mutation rate
const EXT_DOUBLE base_u = global_m/(alph_card-1); // A -> T mutation rate

const char delim = ',';

/* Command-line options: */
uint64_t N = 1E5;
uint64_t T = 1;
uint64_t C = 0;
uint64_t s = 10;

uint64_t verbosity_level = 0;

bool use_fitness_sampler = false;
bool use_simplex_space = false;

bool no_headers = false;
bool skip_writing = false;
bool omit_empties = false;
bool help_flag = false;
bool load_initial_r_from_file = false;
bool randomise_leave_out = true;
bool MLE_exists = false;

NORMAL_DOUBLE GAMMA = -1;
NORMAL_DOUBLE B = 1E-12;

std::string inputFile;

// Global variables:
uint64_t DIM;
uint64_t L;
VectorID counters;
VectorID acc;

uint32_t leave_id;

boost::mutex mutex;

VectorID Data;

MatrixED QT;
MatrixED sQTinv;
MatrixED fQTinv;

MatrixED matrix_f;
MatrixED matrix_m;
MatrixED matrix_s;
MatrixED matrix_p;
MatrixED matrix_r;
MatrixED matrix_ln;

MatrixED population_f;
MatrixED population_m;
MatrixED population_s;
MatrixED population_p;
MatrixED population_r;
MatrixED population_r_new;
VectorED population_ln;
VectorED population_ln_new;

// MCMC options
uint64_t NChains_per_thread;
uint64_t no_samples;

// I/O options
std::string inputFile_initial;

std::string output_f_File("Output-f.txt");

std::string output_m_File("Output-m.txt");
std::string output_s_File("Output-s.txt");

std::string output_p_File("Output-p.txt");
std::string output_r_File("Output-r.txt");

std::string output_ln_File("Output-ln.txt");


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
uint64_t hamming_distance(const std::string& A, const std::string& B)
{
    uint64_t distance = 0;

    for (uint64_t i = 0; i < A.length(); ++i)
    {
        distance += (A[i] != B[i]);
    }

    return distance;
}

void console(const char* str)
{
    std::cout << str << '\n';
}

void console(std::string& str)
{
    console(str.c_str());
}

void about(const char* program_name)
{
    std::cout << "Usage: " << program_name << " [OPTION]... [input-file]\n";
    std::cout << "Infer equilibrium fitness from (presumed) steady-state NGS data\n";
    std::cout << "Example: " << program_name << " -N 1000 -g 0.5 -b 0.1 input.fasta\n\n";


    std::cout << "General options:\n";
    std::cout << "-t            Number of threads to use for sampling.\n";
    std::cout << "              If not specified, defaults to the number of logical cores of the system.\n\n";

    std::cout << "-h, --help    Print this help text\n\n";

    std::cout << "-v[v[v]]      Verbosity level.\n\n";


    std::cout << "MCMC Options:\n";
    std::cout << "-N            Number of total MCMC samples.\n";
    std::cout << "              If not specified, defaults to 10000.\n\n";

    std::cout << "-C            Number of MCMC chains to run in parallel.\n";
    std::cout << "              If not specified, defaults to approximately twice the number of genotypes\n";
    std::cout << "              rounded to a multiple of the number of threads.\n";
    std::cout << "              This ensures an even processing load distribution of the threads.\n\n";

    std::cout << "-s            Chain thinning number.\n";
    std::cout << "              If not specified, defaults to 100\n";
    std::cout << "              Given an MCMC chain, only every s'th element is retained.\n";
    std::cout << "              In general, thinning is required due to the autocorrelation inherent in MCMC.\n\n";

    std::cout << "-g            Differential evolution mixing parameter.\n";
    std::cout << "              If not specified, defaults to FORMULA\n";
    std::cout << "              Differential evolution MCMC requires a scalar for the mixing of the difference vector.\n";
    std::cout << "              This parameter is responsible for achieving optimal target correlation.\n\n";

    std::cout << "-b            Differential evolution jitter parameter.\n";
    std::cout << "              If not specified, defaults to FORMULA\n";
    std::cout << "              Differential evolution MCMC furthermore requires a slight Gaussian jitter in order to\n";
    std::cout << "              guarantee detailed balance. This parameter essentially explores the 'non-Gaussian' part.\n";
    std::cout << "              of the target distribution.\n\n";

    std::cout << "\nI/O Options:\n";
    std::cout << "-o            Omit writing output.\n";
    std::cout << "              This is useful for initial tuning, as the MCMC samples can be quite large.\n\n";

    std::cout << "-0            Omit writing to output genotypes that have no observations.\n";
    std::cout << "              If you have a dataset with 50 observed and 50 non-observed genotypes, this flag will only.\n";
    std::cout << "              print the samples of genotypes with counts larger than 0.\n";

    std::cout << "\n";

    std::cout << "Report bugs to: " << PACKAGE_BUGREPORT << "\n";
    std::cout << PACKAGE_STRING << " home page: <" << PACKAGE_URL << ">\n";
}

void load_inputFile(seq_cont& sequences, seq_inds& indices, VectorED& p_MLE)
{
    typedef std::map<std::string, stats> seq_map;
    typedef std::pair<const std::string, stats> seq_map_pair;

    seq_map all_sequences;

    uint64_t sample_total = 0;

    std::ifstream input;
    input.open(inputFile.c_str());

    if (input.is_open())
    {
        std::string temp;
        stats temp_stats;

        getline(input, temp);
        temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

        input.clear();
        input.seekg(0, std::ios::beg);

        if(temp[0] == '>')
        {
            // FASTA

            if(temp.substr(0, 7) == ">read0_")
            {
                // QuasiRecomb FASTA
                console("Reading QuasiRecomb input file");

                while (input.good())
                {
                    getline(input, temp);
                    temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

                    if (temp.empty())
                        continue;

                    if(temp[0] == '>')
                    {
                        // identifier
                        temp_stats.freq = strtod(temp.substr(temp.find('_')+1, 1000).c_str(), NULL);
                        temp_stats.count = temp_stats.freq*10000;
                        sample_total += temp_stats.count;
                    }
                    else
                    {
                        // sequence
                        all_sequences.insert(seq_map_pair(temp, temp_stats));
                    }
                }
            }
            else
            {
                // Generic FASTA
                console("Reading generic FASTA input file");

                std::string temp_dna;
                std::pair<seq_map::iterator, bool> ret;

                while (input.good())
                {
                    getline(input, temp);
                    temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

                    if (temp.empty())
                        continue;

                    if(temp[0] == '>')
                    {
                        // identifier
                        if (!(temp_dna.empty()))
                        {
                            ret = all_sequences.insert(seq_map_pair(temp_dna, stats()));

                            if (ret.second==false)
                            {
                                ++ret.first->second.count;
                                ++sample_total;
                            }

                            temp_dna.clear();
                        }
                    }
                    else
                    {
                        // sequence
                        temp_dna.append(temp);
                    }
                }

                // add the last line
                if (!(temp_dna.empty()))
                {
                    ret = all_sequences.insert(seq_map_pair(temp_dna, stats()));

                    if (ret.second==false)
                    {
                        ++ret.first->second.count;
                        ++sample_total;
                    }
                }

                for (seq_map::iterator i = all_sequences.begin(); i != all_sequences.end(); ++i)
                {
                    i->second.freq = static_cast<NORMAL_DOUBLE>(i->second.count) / sample_total;
                }
            }
        }
        else
        {
            // QuasiFit input file
            console("Reading generic QuasiFit input file");
            std::string temp_count;

            while (input.good())
            {
                getline(input, temp);
                temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

                if (temp.empty())
                    continue;

                temp_count = temp.substr(temp.find(',')+1,1000);

                if ((omit_empties) && (temp_count == "0"))
                {
                    temp_stats.write_out = false;
                }
                else
                {
                    temp_stats.write_out = true;
                }

                if (temp_count == "K")
                    temp_count = '0';

                temp_stats.count = atoi(temp_count.c_str());
                sample_total += temp_stats.count;

                temp = temp.erase(temp.find(','));

                if (!(temp_stats.count))
                {
                    //need_alpha_pseudocount = true;
                }

                all_sequences.insert(seq_map_pair(temp, temp_stats));
            }

            for (seq_map::iterator i = all_sequences.begin(); i != all_sequences.end(); ++i)
            {
                i->second.freq = static_cast<NORMAL_DOUBLE>(i->second.count) / sample_total;
            }
        }

        input.close();

        DIM = all_sequences.size();

        p_MLE.resize(DIM);
        Data.resize(DIM);

        sequences.reserve(DIM);
        indices.reserve(DIM);

        uint32_t I = 0;
        for (seq_map::iterator i = all_sequences.begin(); i != all_sequences.end(); ++i, ++I)
        {
            sequences.push_back(i->first);

            if (i->second.write_out)
            {
                indices.push_back(I);
            }

            p_MLE(I) = i->second.freq;
            Data(I) = i->second.count;
        }
    }
    else
    {
        console("Input file does not exist");
        exit(EXIT_FAILURE);
    }
}


inline void loadBar()
{
#define BAR_WIDTH 100
#define UPDATE_INTERVAL 1000

    // Only update r times.
    //if ( counter % (NThread/UPDATE_INTERVAL) != 0 )
    //	return;

    NORMAL_DOUBLE ratio = static_cast<NORMAL_DOUBLE>(counters.sum())/(C*N);
    int64_t c = ratio*BAR_WIDTH;
    //std::cout << counter << '\n';

    printf("%5.1f%% [", ratio*100);

    for (int64_t i=0; i<c-1; i++)
        printf("=");
    printf(">");

    for (uint64_t i=c; i<BAR_WIDTH; i++)
        printf(" ");

    printf("]\n\033[F\033[J");
}


uint64_t random_seed()
{
    uint64_t random_seed;
    std::ifstream file ("/dev/urandom", std::ios::binary);
    if (file.is_open())
    {
        char* memblock;
        uint64_t size = sizeof(uint64_t);
        memblock = new char[size];
        file.read(memblock, size);
        file.close();
        random_seed = *reinterpret_cast<uint64_t*>(memblock);
        delete[] memblock;
    }
    else
    {
        std::cout << "Could not open /dev/urandom for generating a seed\n";
        exit(EXIT_FAILURE);
    }
    return random_seed;
}


void counter_display()
{
    while (counters.sum() < N*C)
    {
        loadBar();
        usleep(100E3);
    }
}


void generate_sample(gsl_rng* r, VectorED& new_sample, const uint64_t chain_id)
{
    uint64_t R1, R2;

    // Generate R1
    do
    {
        R1 = gsl_rng_uniform_int(r, C);
    }
    while (R1 == chain_id);

    // Generate R2
    do
    {
        R2 = gsl_rng_uniform_int(r, C);
    }
    while ((R2 == chain_id) || (R2 == R1));

    new_sample = population_r.col(chain_id) + GAMMA*(population_r.col(R1) - population_r.col(R2));

    // Add Gaussian noise
    for (uint32_t i = 0; i < DIM; ++i)
    {
        if (i != leave_id)
            new_sample(i) += gsl_ran_gaussian /*_ziggurat*/ (r, B);
    }
}

EXT_DOUBLE convert_from_probability_to_fitness__manifold(const VectorED& P_vector, VectorED& F_vector, VectorED& M_vector, bool& reject)
{
    EXT_DOUBLE logPosterior = 0;
    M_vector.noalias() = P_vector.asDiagonal().inverse() * fQTinv * P_vector;

    if (M_vector.minCoeff() <= 0)
    {
        reject = true;
    }
    else
    {
        reject = false;

        F_vector.noalias() = M_vector / M_vector.sum();

        Solver LU;
        MatrixED tempM = fQTinv;
        VectorED TempV2 = fQTinv * P_vector;
        MatrixED tempM2 = P_vector.asDiagonal().inverse();
        tempM2 *= TempV2.asDiagonal();
        tempM -= tempM2;

        LU.compute(tempM);

        for (uint64_t i = 0; i < DIM; ++i)
        {
            logPosterior += log(fabs(LU.matrixLU()(i, i))) + Data(i)*log(P_vector(i));
        }
    }

    return logPosterior;
}

inline void convert_from_M_to_S(const VectorED& M_vector, VectorED& S_vector)
{
    S_vector = M_vector - VectorED::Ones(DIM);
}

EXT_DOUBLE (*fitness_space) (const VectorED& P_vector, VectorED& F_vector, VectorED& M_vector, bool& reject) = convert_from_probability_to_fitness__manifold;

void convert_from_R_to_P(const VectorED& R_vector, VectorED& P_vector_unnorm, VectorED& P_vector)
{
    EXT_DOUBLE sum = 0;
    //P_vector_unnorm(DIM-1) = 1;

    for (uint32_t i = 0; i < DIM; ++i)
    {
        if (i == leave_id)
            P_vector_unnorm(i) = 1;
        else
            P_vector_unnorm(i) = exp(R_vector(i));

        sum += P_vector_unnorm(i);
    }

    P_vector = P_vector_unnorm / sum;
}

void convert_from_P_to_R(const VectorED& P_vector, VectorED& R_vector)
{
    for (uint32_t i = 0; i < DIM; ++i)
    {
        if (i == leave_id)
            R_vector(i) = 0;
        else
            R_vector(i) = log(P_vector(i) / P_vector(leave_id));
    }
}

// 2. probability sampler
void probability_sampler(
    boost::barrier& syn_barrier,
    uint64_t thread_no)
{
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    //gsl_rng* r = gsl_rng_alloc(gsl_rng_ranlxd2);
    //gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
    uint64_t seed;
    seed = random_seed();
    gsl_rng_set(r, seed);

    bool reject_outright;

    if (thread_no == 0)
    {
        std::cout << "Probability sampler\n";
        std::cout << "gamma: " << GAMMA << '\n';
        std::cout << "b:     " << B << '\n';
        if (!(randomise_leave_out))
            std::cout << "NOT randomising position of leave-one-out in r vector\n";
    }

    EXT_DOUBLE temp_logPosterior;

    VectorED r_new(DIM);
    VectorED p_new(DIM), p_new_unnorm(DIM), m_new(DIM), s_new(DIM), f_new(DIM);

    for (uint64_t i = 0; i < N; ++i)
    {
        // Perform N MCMC trials ...

        for (uint64_t j = thread_no*NChains_per_thread; j < (thread_no+1)*NChains_per_thread; ++j)
        {
            // ... for NChains_per_thread Chains
            ++counters(j);
            //reject_outright = false;

            // 1) generate new R
            generate_sample(r, r_new, j);

            // 2) convert R to P + unnormalised P
            convert_from_R_to_P(r_new, p_new_unnorm, p_new);

            // 3) convert P to M and F
            temp_logPosterior = fitness_space(p_new, f_new, m_new, reject_outright);

            if (verbosity_level >= 4)
            {
                mutex.lock();
                std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
                std::cout << "r_new: " << r_new.transpose().cast<long double>() << '\n';
                std::cout << "m_new: " << m_new.transpose().cast<long double>() << '\n';
                std::cout << "p_new: " << p_new.transpose().cast<long double>() << '\n';
                mutex.unlock();
            }

            if ((!(reject_outright)) && ((temp_logPosterior - population_ln(0, j)) > log(gsl_ran_flat(r, 0, 1))))
            {
                // accept proposal
                // 4) convert M to S
                convert_from_M_to_S(m_new, s_new);

                population_ln_new(0, j) = temp_logPosterior;
                population_r_new.col(j) = r_new;

                population_p.col(j) = p_new;
                population_f.col(j) = f_new;
                population_s.col(j) = s_new;
                population_m.col(j) = m_new;

                ++acc(j);

                if ((verbosity_level >= 3) && (m_new(1) < 1))
                {
                    mutex.lock();
                    std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
                    std::cout << "r_new: " << r_new.transpose().cast<long double>() << '\n';
                    std::cout << "m_new: " << m_new.transpose().cast<long double>() << '\n';
                    std::cout << "p_new: " << p_new.transpose().cast<long double>() << '\n';
                    mutex.unlock();
                }
            }
            else
            {
                // reject proposal
                population_ln_new(0, j) = population_ln(0, j);
                population_r_new.col(j) = population_r.col(j);
            }
        }
        syn_barrier.wait();

        if (thread_no == 0)
        {
            population_r.swap(population_r_new);
            population_ln.swap(population_ln_new);

            if (i % s == 0)
            {
                matrix_f.block(0, i/s*C, DIM, C) = population_f;
                matrix_m.block(0, i/s*C, DIM, C) = population_m;
                matrix_s.block(0, i/s*C, DIM, C) = population_s;
                matrix_p.block(0, i/s*C, DIM, C) = population_p;
                matrix_r.block(0, i/s*C, DIM, C) = population_r;
                matrix_ln.block(0, i/s*C, 1, C) = population_ln;
            }

            if (randomise_leave_out)
                leave_id = gsl_rng_uniform_int(r, DIM);
        }
        syn_barrier.wait();

        // regenerate r population
        if (randomise_leave_out)
        {
            for (uint64_t j = thread_no*NChains_per_thread; j < (thread_no+1)*NChains_per_thread; ++j)
            {
                convert_from_P_to_R(population_p.col(j), r_new);
                population_r.col(j) = r_new;
            }
            syn_barrier.wait();
        }
    }

    gsl_rng_free(r);
}

void (*MCMC)(
    boost::barrier& syn_barrier,
    uint64_t thread_no) = probability_sampler;


void r_loader_from_file()
{
    std::ifstream input;
    std::string temp;

    std::vector<std::string> line;

    input.open(inputFile_initial.c_str());
    uint32_t i = 0;

    if (input.is_open())
    {
        std::cout << "Loading initial population from " << inputFile_initial << '\n';

        while (input.good())
        {
            getline(input, temp);
            temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());
            if (temp.empty())
                continue;

            boost::split(line, temp, boost::is_any_of(","), boost::token_compress_on);

            if (line.size() != DIM)
            {
                std::cout << "FATAL ERROR: Initial number of dimensions do not match!\n";
                std::cout << "Expected number of comma separated items: " << DIM << '\n';
                std::cout << "  Actual number of comma separated items: " << line.size() << '\n';
                exit(EXIT_FAILURE);
            }

            for (uint32_t j = 0; j < DIM; ++j)
            {
                population_r(j, i) = strtod(line[j].c_str(), NULL);
            }

            ++i;
        }


        std::cout << "Read " << i << " lines from " << inputFile_initial;

        if (i != C)
        {
            std::cout << ", but expected " << C << '\n';
            exit(EXIT_FAILURE);
        }

    }
    else
    {
        std::cout << inputFile_initial << " does not exist, aborting\n";
        exit(EXIT_FAILURE);
    }

    for (uint32_t i = 0; i < DIM; ++i)
    {
        if (population_r.row(i) == VectorED::Zero(C))
            leave_id = i;
    }

    std::cout << "\nDetected leave out column to be " << leave_id << '\n';
    input.close();
}


void population_initialiser(uint64_t thread_no)
{
    gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
    uint64_t seed;
    seed = random_seed();
    gsl_rng_set(r, seed);

    VectorED r_new(DIM);
    VectorED f_new(DIM);
    VectorED p_new(DIM), p_new_unnorm(DIM);
    VectorED s_new(DIM), m_new(DIM);
    EXT_DOUBLE temp_logPosterior;
    bool reject_outright;

    for (int j = C*thread_no/T; j < C*(thread_no+1)/T; ++j)
    {
        // 1) convert R to P + unnormalised P
        convert_from_R_to_P(population_r.col(j), p_new_unnorm, p_new);

        // 2) convert P to M and F
        temp_logPosterior = fitness_space(p_new, f_new, m_new, reject_outright);

        // 3) convert M to S
        convert_from_M_to_S(m_new, s_new);

        population_ln(0, j) = temp_logPosterior;

        population_p.col(j) = p_new;
        population_f.col(j) = f_new;
        population_s.col(j) = s_new;
        population_m.col(j) = m_new;
    }

    gsl_rng_free(r);
}


void parse_arguments(int argc, char** argv)
{
    int c, option_index = 0;

    while ((c = getopt_long (argc, argv, "t:N:C:s:g:b:Ro0vh", long_options, &option_index)) != -1)
    {
        switch (c)
        {
        case 't':
            T = atoi(optarg);
            break;

            // MCMC options
        case 'N':
            N = strtoull(optarg, NULL, 10);
            break;

        case 'C':
            C = strtoull(optarg, NULL, 10);
            break;

        case 's':
            s = atoi(optarg);
            break;

        case 1000:
            use_fitness_sampler = true;
            break;

        case 1001:
            use_simplex_space = true;
            break;

        case 'g':
            GAMMA = strtod(optarg, NULL);
            break;

        case 'b':
            B = strtod(optarg, NULL);
            break;

        case 1020:
            randomise_leave_out = false;
            break;

        case 1030:
            load_initial_r_from_file = true;
            break;

            // Input/Output options
        case 1010:
            no_headers = true;
            break;

        case 'o':
            skip_writing = true;
            break;

        case '0':
            omit_empties = true;
            break;

            // Miscellaneous options
        case 'v':
            ++verbosity_level;
            break;

        case 'h':
            help_flag = true;
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf ("option %s", long_options[option_index].name);
            if (optarg)
                printf (" with arg %s", optarg);
            printf ("\n");
            break;

        default:
            abort ();
        }
    }

    if ((help_flag) || (argc == 1))
    {
        about(argv[0]);
        exit(EXIT_SUCCESS);
    }

    if (optind < argc)
    {
        inputFile = argv[optind];
    }
    else
    {
        std::cout << "Missing input file\n";
        exit(EXIT_FAILURE);
    }

    if (use_fitness_sampler && use_simplex_space)
    {
        std::cout << "Either use the fitness sampler or employ the probability sampler with fitness space in the simplex\n";
        exit(EXIT_FAILURE);
    }
}
