#include <quasifit.hpp>

/* Custom console */
console_buf console_buffer;
std::ostream console(&console_buffer);

/* Input handling variables */
const char delim = ',';
std::string inputFile;

/* Command-line options */
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

/* Global variables */
uint64_t DIM;
uint64_t Nreads;
uint64_t L;

VectorID Data;

// I/O options
std::string inputFile_initial;
std::string output_f_File("Output-f.txt");
std::string output_m_File("Output-m.txt");
std::string output_s_File("Output-s.txt");
std::string output_p_File("Output-p.txt");
std::string output_r_File("Output-r.txt");
std::string output_ln_File("Output-ln.txt");


uint64_t hamming_distance(const std::string& A, const std::string& B)
{
    uint64_t distance = 0;

    for (uint64_t i = 0; i < A.length(); ++i)
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
    std::cout << PACKAGE_STRING << " built with " << PRECISION_TYPE << " (~" << std::numeric_limits<EXT_DOUBLE>::digits10 << " digits)\n";
    std::cout << "Home page: <" << PACKAGE_URL << ">\n";
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
                console << "Reading QuasiRecomb input file\n";

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
                console << "Reading generic FASTA input file\n";

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
            console << "Reading generic QuasiFit input file\n";
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

            Data(I) = i->second.count;
        }

        Nreads = Data.sum();
        if (Nreads)
        {
            p_MLE = Data.cast<EXT_DOUBLE>() / Nreads;
        }
        else
        {
            p_MLE = VectorED::Constant(DIM, 1.0/DIM);
        }

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
    ez::ezETAProgressBar p(C*N-1);
    p.start();

    do
    {
        current = counters.sum();
        p = current;
        usleep(100E3);
    }
    while (counters.sum() < N*C);
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
