/* Command-line options */
extern uint32_t verbosity_level;

extern bool no_headers;
extern bool skip_writing;
extern bool help_flag;
extern bool load_initial_r_from_file;
extern bool randomise_leave_out;



/* Data-containing variables */
extern VectorID Data;
extern VectorID closest_observed_neighbour;

extern VectorED p_MLE;
extern bool MLE_exists;

extern uint32_t DIM;
extern uint32_t Nreads;
extern uint32_t L;



/* MCMC parameters */
extern uint64_t N;
extern int32_t T;
extern int32_t C;
extern uint32_t s;

extern uint32_t NChains_per_thread;
extern uint64_t no_samples;

extern NORMAL_DOUBLE GAMMA;
extern NORMAL_DOUBLE B;



/* Thread variables */
extern uint32_t leave_id;
extern boost::mutex mutex;
extern VectorID counters;
extern VectorID shuffled_index;
extern VectorID acc;



/* Matrices/Vectors/Constants for calculating transformations */
extern const uint32_t alph_card;
extern const EXT_DOUBLE global_m;	// global HIV mutation rate
extern const EXT_DOUBLE base_u;	// A -> T mutation rate

extern MatrixED QT;
/*
   extern MatrixED sQTinv;
   extern MatrixED fQTinv;
 */

extern Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;
extern MatrixED Q_Random;



/* Storage variables */
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



/* I/O variables */
extern const char delim;

extern std::ostream console;

extern std::string inputFile;

extern std::string inputFile_initial;
extern std::string output_f_File;
extern std::string output_m_File;
extern std::string output_s_File;
extern std::string output_p_File;
extern std::string output_r_File;
extern std::string output_ln_File;