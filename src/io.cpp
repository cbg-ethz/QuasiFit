#include <quasifit.hpp>
#include <getopt.h>

/* Custom console */
console_buf console_buffer;
std::ostream console(&console_buffer);

/* Input handling variables */
const char delim = ',';
std::string inputFile;

/* Command-line options */
uint64_t N = 100000;
int32_t T = -1;
int32_t C = -1;
uint32_t s = 25;

uint32_t verbosity_level = 0;

bool no_headers = false;
bool skip_writing = false;
bool help_flag = false;
bool load_initial_r_from_file = false;
bool randomise_leave_out = true;
bool MLE_exists = true;

NORMAL_DOUBLE GAMMA = -1;
NORMAL_DOUBLE B = 1E-12;

/* Global variables */
uint32_t DIM;
uint32_t Nreads;
uint32_t L;

VectorID Data;
VectorID closest_observed_neighbour;

// I/O options
std::string inputFile_initial;
std::string output_f_File;
std::string output_p_File;
std::string output_r_File;
std::string output_ln_File;


uint32_t hamming_distance(const std::string& A, const std::string& B)
{
	uint32_t distance = 0;

	for (uint32_t i = 0; i < A.length(); ++i)
	{
		distance += (A[i] != B[i]);
	}

	return distance;
}

void about(const char* program_name)
{
	std::cout << "Usage: " << program_name << " [OPTION]... [input-file]\n";
	std::cout << "Infer equilibrium fitness from (presumed) steady-state NGS data\n";
	std::cout << "Example: " << program_name << " -N 1000 -g 0.5 -b 0.1 input.fasta\n\n";

	std::cout << "General options:\n";
	std::cout << "-t            Number of threads to use for sampling\n";
	std::cout << "              If not specified, defaults to the number of logical cores of the system.\n\n";
	std::cout << "-h, --help    Print this help text\n\n";
	std::cout << "-v[v[v]]      Verbosity level\n\n";

	std::cout << "MCMC Options:\n";
	std::cout << "-N            Number of total MCMC samples PER CHAIN\n";
	std::cout << "              If not specified, defaults to 100'000.\n\n";
	std::cout << "-C            Number of MCMC chains to run in parallel\n";
	std::cout << "              If not specified, defaults to approximately twice the number of genotypes\n";
	std::cout << "              rounded to a multiple of the number of threads.\n";
	std::cout << "              This ensures an even processing load distribution of the threads across cores.\n\n";
	std::cout << "-s            Chain thinning number\n";
	std::cout << "              If not specified, defaults to 25.\n";
	std::cout << "              Given an MCMC chain, only every s'th element is retained.\n";
	std::cout << "              In general, thinning is required due to the autocorrelation inherent in MCMC.\n\n";
	std::cout << "-g            Differential evolution mixing parameter\n";
	std::cout << "              If not specified, defaults to min(1.6829 / sqrt(n), 0.7), where n is the dimensionality.\n";
	std::cout << "              Differential evolution MCMC requires a scalar for the mixing of the difference vector.\n";
	std::cout << "              This parameter is responsible for achieving optimal target correlation.\n\n";
	std::cout << "-b            Differential evolution jitter parameter\n";
	std::cout << "              If not specified, defaults to 1E-12.\n";
	std::cout << "              Differential evolution MCMC furthermore requires a slight Gaussian jitter in order to\n";
	std::cout << "              guarantee detailed balance. This parameter essentially explores the 'non-Gaussian' part\n";
	std::cout << "              of the target distribution.\n\n";
	std::cout << "--nR          Do not randomise the index to be used as normalising factor in the transformation t(y).\n";
	std::cout << "              It is STRONGLY advised against setting this option, as numerical errors will not be\n";
	std::cout << "              distributed randomly.\n\n";

	std::cout << "I/O Options:\n";
	std::cout << "-o            Omit writing output\n";
	std::cout << "              This is useful for initial tuning, as the MCMC samples can be quite large.\n\n";
	std::cout << "--nH          Do not write the header line for output files\n";
	std::cout << "              This is useful if you want to load the output files into MATLAB.\n\n";

	std::cout << "Report bugs to: " << PACKAGE_BUGREPORT << "\n";
	std::cout << PACKAGE_STRING << " built with " << PRECISION_TYPE << " (~" << std::numeric_limits<EXT_DOUBLE>::digits10 << " digits)\n";
	std::cout << "Home page: <" << PACKAGE_URL << ">\n";
}

void load_inputFile(seq_cont& sequences)
{
	typedef std::map<std::string, stats> seq_map;
	typedef std::pair<const std::string, stats> seq_map_pair;

	// input_type = 0: Generic FASTA
	// input_type = 1: QuasiRecomb FASTA
	// input_type = 2: QuasiFit
	char input_type;

	seq_map all_sequences;

	std::ifstream input;
	input.open(inputFile.c_str());

	if (input.is_open())
	{
		std::string temp;
		stats temp_stats;

		getline(input, temp);
		temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

		// 1. First detect type of input
		if (temp[0] == '>')	// FASTA
		{
			if (temp.substr(0, 7) == ">read0_")
			{
				// QuasiRecomb FASTA
				console << "Reading QuasiRecomb input file\n";
				input_type = 1;
			}
			else
			{
				// Generic FASTA
				console << "Reading generic FASTA input file\n";
				input_type = 0;
			}
		}
		else
		{
			// QuasiFit input file
			console << "Reading generic QuasiFit input file\n";
			input_type = 2;
		}

		input.clear();
		input.seekg(0, std::ios::beg);

		// 2. Now read specific type
		switch (input_type)
		{
			case 0:	// Generic FASTA
			{
				std::string temp_dna;
				std::pair<seq_map::iterator, bool> ret;

				while (input.good())
				{
					getline(input, temp);
					temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

					if (temp.empty())
						continue;

					if (temp[0] == '>')
					{
						// identifier
						if (!(temp_dna.empty()))
						{
							ret = all_sequences.insert(seq_map_pair(temp_dna, stats()));

							if (ret.second == false)
								++ret.first->second.count;

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

					if (ret.second == false)
					{
						++ret.first->second.count;
					}
				}
				break;
			}

			case 1:	// QuasiRecomb FASTA
			{
				while (input.good())
				{
					getline(input, temp);
					temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

					if (temp.empty())
						continue;

					if (temp[0] == '>')
					{
						// identifier
						temp_stats.count = strtod(temp.substr(temp.find('_') + 1, 1000).c_str(), NULL) * 10000;
					}
					else
					{
						// sequence
						all_sequences.insert(seq_map_pair(temp, temp_stats));
					}
				}
				break;
			}

			case 2:	// QuasiFit
			{
				std::string temp_count;

				while (input.good())
				{
					getline(input, temp);
					temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());

					if (temp.empty())
						continue;

					temp_count = temp.substr(temp.find(',') + 1,1000);
					temp_stats.count = atoi(temp_count.c_str());

					temp = temp.erase(temp.find(','));
					all_sequences.insert(seq_map_pair(temp, temp_stats));
				}
				break;
			}
		}

		input.close();

		DIM = all_sequences.size();

		Data.resize(DIM);
		sequences.reserve(DIM);

		uint32_t I = 0;
		L = all_sequences.begin()->first.length();
		for (seq_map::iterator i = all_sequences.begin(); i != all_sequences.end(); ++i, ++I)
		{
			sequences.push_back(i->first);
			Data(I) = i->second.count;

			if (L != i->first.length())
			{
				console << "Haplotype " << i->first << " has a different length than the rest, aborting\n";
				exit(EXIT_FAILURE);
			}
		}

		Nreads = Data.sum();
	}
	else
	{
		console << "Input file does not exist\n";
		exit(EXIT_FAILURE);
	}
}

void counter_display()
{
	uint64_t current;
	ez::ezETAProgressBar p(C * N - 1);
	p.start();

	do
	{
		current = counters.sum();
		p = current;
		usleep(100E3);
	}
	while (counters.sum() < N * C);
}

void r_loader_from_file()
{
	std::ifstream input;
	std::string temp;

	std::vector<std::string> line;

	input.open(inputFile_initial.c_str());
	uint32_t i = 0;

	if (input.is_open())
	{
		console << "Loading initial population from " << inputFile_initial << '\n';

		while (input.good())
		{
			getline(input, temp);
			temp.erase(std::remove(temp.begin(), temp.end(), '\n'), temp.end());
			if (temp.empty())
				continue;

			boost::split(line, temp, boost::is_any_of(","), boost::token_compress_on);

			if (line.size() != DIM)
			{
				console << "FATAL ERROR: Initial number of dimensions do not match!\n";
				console << "Expected number of comma separated items: " << DIM << '\n';
				console << "  Actual number of comma separated items: " << line.size() << '\n';
				exit(EXIT_FAILURE);
			}

			for (uint32_t j = 0; j < DIM; ++j)
			{
				population_r(j, i) = strtod(line[j].c_str(), NULL);
			}

			++i;
		}


		console << "Read " << i << " lines from " << inputFile_initial;

		if (i != C)
		{
			console << ", but expected " << C << '\n';
			exit(EXIT_FAILURE);
		}

	}
	else
	{
		console << inputFile_initial << " does not exist, aborting\n";
		exit(EXIT_FAILURE);
	}

	for (uint32_t i = 0; i < DIM; ++i)
	{
		if (population_r.row(i) == VectorED::Zero(C))
			leave_id = i;
	}

	console << "\nDetected leave out column to be " << leave_id << '\n';
	input.close();
}

static struct option long_options[] =
{
	{"nH", no_argument,       0, 1010},
	{"nR", no_argument,       0, 1020},
	{"l",  required_argument, 0, 1030},
	{0, 0, 0, 0}
};

void parse_arguments(int argc, char** argv)
{
	int c, option_index = 0;

	while ((c = getopt_long (argc, argv, "t:N:C:s:g:b:ovh", long_options, &option_index)) != -1)
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
}

void write_output(const seq_cont& sequences)
{
	/* Output the samples */
	std::ofstream output;

	if (skip_writing)
	{
		console << "Skipped writing output files\n";
	}
	else
	{
		Eigen::IOFormat QuasiFitFormat_F(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n", "", "", "", "");
		Eigen::IOFormat QuasiFitFormat(8, Eigen::DontAlignCols, ",", "\n", "", "", "", "");

		if (no_headers)
			console << "Writing no headers to output files\n";

		/* fitness manifold samples */
		output.open(output_f_File.c_str());
		// header:
		if (!(no_headers))
		{
			output << sequences[0];
			for (uint32_t j = 1; j < DIM; ++j)
			{
				output << delim << sequences[j];
			}
			output << '\n';
		}
		// actual data:
		output << matrix_f.transpose().format(QuasiFitFormat_F);
		output.close();


		/* probability samples */
		output.open(output_p_File.c_str());
		// header:
		if (!(no_headers))
		{
			output << sequences[0];
			for (uint32_t j = 1; j < DIM; ++j)
			{
				output << delim << sequences[j];
			}
			output << '\n';
		}
		// actual data:
		output << matrix_p.transpose().format(QuasiFitFormat);
		output.close();


		/* euclidean samples */
		output.open(output_r_File.c_str());
		// header:
		if (!(no_headers))
		{
			output << sequences[0];
			for (uint32_t j = 1; j < DIM; ++j)
			{
				output << delim << sequences[j];
			}
			output << '\n';
		}
		// actual data:
		output << matrix_r.transpose().format(QuasiFitFormat);
		output.close();


		/* diagnostic samples */
		output.open(output_ln_File.c_str());
		// header:
		if (!(no_headers))
		{
			output << "LogPost,LogDet,LogMult\n";
		}
		// actual data:
		output << matrix_ln.transpose().format(QuasiFitFormat);
		output.close();
	}
}