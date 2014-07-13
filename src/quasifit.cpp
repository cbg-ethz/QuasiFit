#include <quasifit.hpp>

/* mutation rates */
EXT_DOUBLE mu = 3E-5;	// global HIV mutation rate
EXT_DOUBLE kappa = 4;	// transition/transversion rate

void initialise_matrices(const seq_cont& seqs);

int main (int argc, char** argv)
{
	/* 1) parse command line and greeting */
	parse_arguments(argc, argv);
	console << "This is " << PACKAGE_STRING << '\n';
	console << "Home page: <" << PACKAGE_URL << ">\n";

	/* 2) process input file */
	seq_cont sequences;
	load_inputFile(sequences);

	/* 3) Setup filenames */
	size_t pos = inputFile.find_last_of('.');
	if (pos != std::string::npos)
	{
		output_f_File = inputFile.substr(0, pos) + "-f.csv";
		output_p_File = inputFile.substr(0, pos) + "-p.csv";
		output_r_File = inputFile.substr(0, pos) + "-r.csv";
		output_ln_File = inputFile.substr(0, pos) + "-diag.csv";
	}
	else
	{
		output_f_File = inputFile + "-f.csv";
		output_p_File = inputFile + "-p.csv";
		output_r_File = inputFile + "-r.csv";
		output_ln_File = inputFile + "-diag.csv";
	}

	/* 4) initialize matrices for h(p) */
	initialise_matrices(sequences);

	/* 5) setup and check parameters */
	if (GAMMA == -1)
		GAMMA = std::min(1.6829 / sqrt(DIM), 0.7);

	if (T == -1)
		T = boost::thread::hardware_concurrency();

	if (C == -1)
		C = 2 * DIM;

	bool reducedThreads = false;
	if (T > C)
	{
		T = C;
		reducedThreads = true;
	}

	NChains_per_thread = ceil(1.0 * C / T);
	C = NChains_per_thread * T;

	no_samples = N / s;
	N = no_samples * s;

	/* 6) allocate memory for matrices holding samples */
	matrix_f.resize(DIM, no_samples * C);
	matrix_p.resize(DIM, no_samples * C);
	matrix_r.resize(DIM, no_samples * C);
	matrix_ln.resize(3, no_samples * C);

	population_f.resize(DIM, C);
	population_p.resize(DIM, C);
	population_r.resize(DIM, C);
	population_r_new.resize(DIM, C);
	population_ln.resize(3, C);
	population_ln_new.resize(3, C);

	counters.resize(C);
	counters.setZero();

	acc.resize(C);
	acc.setZero();

	/* 7) Initialise population */
	population_initialiser();

	/* 8) Print diagnostic information */
	console << "Number of haplotypes:   " << DIM << '\n';
	if (!(Nreads))
		console << "Sampling from prior\n";
	else
		console << "Number of total counts: " << Nreads << '\n';

	if (MLE_exists)
		console << "MLE exists\n";
	else
		console << "MLE does NOT exist - expect slower convergence and decreased efficiency\n";

	if (verbosity_level >= 2)
	{
		for (uint32_t i = 0; i < DIM; ++i)
		{
			console << sequences[i] << '\t' << Data(i) << '\t' << p_MLE(i) << '\n';
		}
	}

	if (reducedThreads)
		console << "Reduced number of threads to:   " << T << '\n';
	else
		console << "Number of threads:              " << T << '\n';

	console << "Number of trials (N) per chain: " << N << '\n';
	console << "Number of trials in total:      " << N * T << '\n';
	console << "Number of chains (C):           " << C << '\n';
	console << "Chains per thread:              " << NChains_per_thread << '\n';
	console << "Thinning number (s):            " << s << '\n';
	console << "Number of samples per chain:    " << no_samples << '\n';
	console << "Number of samples in total:     " << no_samples * C << '\n';

	console << "Gamma: " << GAMMA << '\n';
	console << "B:     " << B << '\n';
	if (!(randomise_leave_out))
	{
		console << "NOT randomising position of leave-one-out in r vector\n";
		console << "Tread carefully, numerical issues abound\n";
	}

	/* 9) Start MCMC threads */
	boost::chrono::process_real_cpu_clock::time_point realEnd, realStart = boost::chrono::process_real_cpu_clock::now();
	boost::chrono::process_user_cpu_clock::time_point userEnd, userStart = boost::chrono::process_user_cpu_clock::now();
	boost::chrono::process_system_cpu_clock::time_point systemEnd, systemStart = boost::chrono::process_system_cpu_clock::now();
	boost::chrono::high_resolution_clock::time_point timeEnd, timeStart = boost::chrono::high_resolution_clock::now();

	boost::thread_group* all_threads = new boost::thread_group;
	boost::barrier syn_barrier(T);

	if (T > 5)
		console << "Starting MCMC threads\n";

	for (uint32_t t = 0; t < T; ++t)
	{
		if (T <= 5)
			console << "Starting thread #" << t << '\n';

		all_threads->add_thread(new boost::thread(MCMC, boost::ref(syn_barrier), t));
	}

	if (verbosity_level >= 1)
		all_threads->add_thread(new boost::thread(counter_display));

	console << "Processing all threads\n";
	all_threads->join_all();
	delete all_threads;

	realEnd = boost::chrono::process_real_cpu_clock::now();
	userEnd = boost::chrono::process_user_cpu_clock::now();
	systemEnd = boost::chrono::process_system_cpu_clock::now();

	timeEnd = boost::chrono::high_resolution_clock::now();

	typedef boost::chrono::duration<int, boost::milli> millisecs_t;
	typedef boost::chrono::duration<double> secs_t;

	secs_t lReal(boost::chrono::duration_cast<millisecs_t>(realEnd - realStart)),
	lUser(boost::chrono::duration_cast<millisecs_t>(userEnd - userStart)),
	lSystem(boost::chrono::duration_cast<millisecs_t>(systemEnd - systemStart)),
	lTime(boost::chrono::duration_cast<millisecs_t>(timeEnd - timeStart));

	console << std::fixed << std::setprecision(3);
	console << "MCMC ran for:\n";
	console << "     Total Time: " << lTime.count() << "s\n";
	console << "  Real CPU Time: " << lReal.count() << "s\n";
	console << "  User CPU Time: " << lUser.count() << "s\n";
	console << "System CPU Time: " << lSystem.count() << "s\n";

	console << "Statistics:\n";
	console << "Acceptance rate: " << static_cast<NORMAL_DOUBLE>(acc.sum()) / (N * C) * 100 << "%\n";

	/* 10) Write output files */
	write_output(sequences);
	return 0;
}

void initialise_matrices(const seq_cont& seqs)
{
	MatrixID H_alpha(DIM, DIM);
	MatrixID H_beta(DIM, DIM);
	H_alpha.setZero();
	H_beta.setZero();

	int n_ti, n_tv;
	unsigned char sum;

	// set up both hamming matrices
	for (uint32_t i = 0; i < DIM - 1; ++i)
	{
		for (uint32_t j = i + 1; j < DIM; ++j)
		{
			n_ti = 0;
			n_tv = 0;

			for (int k = 0; k < L; ++k)
			{
				if (seqs[i][k] != seqs[j][k])
				{
					// mutation
					sum = seqs[i][k] + seqs[j][k];

					if ((sum == 136) || (sum == 151))
					{
						// transition
						++n_ti;
					}
					else
					{
						// transversion
						++n_tv;
					}
				}
			}

			H_alpha(i, j) = n_ti;
			H_alpha(j, i) = n_ti;

			H_beta(i, j) = n_tv;
			H_beta(j, i) = n_tv;
		}
	}

	EXT_DOUBLE alpha = mu * kappa / (kappa + 2);
	EXT_DOUBLE beta = mu / (kappa + 2);
	console << "Mutation rate (mu):  " << mu << '\n';
	console << "ti/tv rate (kappa):  " << kappa << '\n';

	if (kappa == 1)
	{
		console << "Using standard uniform quasispecies mutation model.\n";
	}
	else
	{
		console << "Transition (alpha):  " << alpha << '\n';
		console << "Transversion (beta): " << beta << '\n';
	}

	if (verbosity_level >= 3)
	{
		if (kappa == 1)
			console << "Hamming Matrix:\n" << H_alpha + H_beta << '\n';
		else
			console << "Transition Hamming Matrix:\n" << H_alpha << "\nTransversion Hamming Matrix:\n" << H_beta << '\n';
	}

	// set up mutation matrix Q
	QT = MatrixED::Identity(DIM, DIM);
	closest_observed_neighbour.resize(DIM);

	EXT_DOUBLE q_ij, q_cum, q_max;
	std::vector<EXT_DOUBLE> q_terms;
	for (uint32_t i = 0; i < DIM; ++i)
	{
		q_terms.clear();
		q_max = 0;

		for (uint32_t j = 0; j < DIM; ++j)
		{
			if (i != j)
			{
				q_ij = pow(alpha, H_alpha(i, j)) * pow(beta, H_beta(i, j)) * pow(1 - mu, L - H_alpha(i, j) - H_beta(i, j));
				q_terms.push_back(q_ij);

				// findest closest observed mutation partner
				if ((q_ij > q_max) && (Data(j)))
				{
					q_max = q_ij;
				}

				QT(i, j) = q_ij;
			}
		}

		closest_observed_neighbour(i) = q_max;

		// collect all terms to calculate seld-replication probability
		std::sort(q_terms.begin(), q_terms.end());
		q_cum = 0;

		for (std::vector<EXT_DOUBLE>::const_iterator j = q_terms.begin(); j != q_terms.end(); ++j)
			q_cum += *j;

		QT(i, i) -= q_cum;
	}

	QT.transposeInPlace();

	// Check assertions
	for (uint32_t i = 0; i < DIM; ++i)
	{
		for (uint32_t j = i; j < DIM; ++j)
		{
			// Symmetry
			if (QT(i, j) != QT(j, i))
			{
				console << "Matrix Q is not symmetrical!\n";
				exit(EXIT_FAILURE);
			}
		}
	}

	P.resize(DIM);
	P.setIdentity();
	std::random_shuffle(P.indices().data(), P.indices().data() + P.indices().size());

	Q_Random.noalias() = P * QT * P.transpose();

	// fill shuffling vector
	shuffled_index.resize(DIM);
	for (uint32_t i = 0; i < DIM; ++i)
	{
		shuffled_index[i] = i;
	}

	if (verbosity_level >= 3)
	{
		console << std::scientific << std::setprecision(6);
		console << "QT:\n" << QT << '\n';
	}
}

