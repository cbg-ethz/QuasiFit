#include <quasifit.hpp>

/* Constants */
const uint32_t alph_card = 4;
const EXT_DOUBLE global_m = 3E-5;	// global HIV mutation rate
const EXT_DOUBLE base_u = global_m / (alph_card - 1);	// A -> T mutation rate

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

	if (T > C)
	{
		T = C;
		console << "Reduced number of threads to " << T << '\n';
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
		console << "MLE does NOT exist - expect slower convergence and efficiency\n";

	if (verbosity_level >= 2)
	{
		for (uint32_t i = 0; i < DIM; ++i)
		{
			console << sequences[i] << '\t' << Data(i) << '\t' << p_MLE(i) << '\n';
		}
	}

	console << "Number of threads:              " << T << '\n';
	console << "Number of trials (N) per chain: " << N << '\n';
	console << "Number of chains (C):           " << C << '\n';
	console << "Chains per thread:              " << NChains_per_thread << '\n';
	console << "Thinning number (s):            " << s << '\n';
	console << "Number of samples in total:     " << no_samples << '\n';

	console << "Gamma: " << GAMMA << '\n';
	console << "B:     " << B << '\n';
	if (!(randomise_leave_out))
	{
		console << "NOT randomising position of leave-one-out in r vector\n";
		console << "Tread carefully, numerical issues abound\n";
	}

	/* 9) Start MCMC threads*/
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

	write_output(sequences);
	return 0;
}

void initialise_matrices(const seq_cont& seqs) {
	MatrixID H(DIM, DIM);
	std::vector<VectorID> MutationNeighbours(DIM);
	VectorID temp_Profile(L);
	closest_observed_neighbour.resize(DIM);
	int32_t dist;

	// create Hamming distance matrix
	for (uint32_t i = 0; i < DIM; ++i)
	{
		temp_Profile.setZero();
		closest_observed_neighbour(i) = L;

		for (uint32_t j = 0; j < DIM; ++j)
		{
			dist = hamming_distance(seqs[i], seqs[j]);
			H(i, j) = dist;
			if (i != j)
			{
				++temp_Profile(H(i, j) - 1);
				if ((Data(j)) && (dist < closest_observed_neighbour(i)))
					closest_observed_neighbour(i) = dist;
			}
		}

		MutationNeighbours[i] = temp_Profile;
	}

	if (verbosity_level >= 3)
		console << H << '\n';

	// Setup Q Matrix
	MatrixED D(DIM, DIM);
	EXT_DOUBLE q_cum, q_ij;

	for (uint64_t i = 0; i < DIM; ++i)
	{
		// fill Q matrix generically
		for (uint64_t j = 0; j < DIM; ++j)
		{
			if (i != j)
			{
				q_ij = pow(base_u, H(i, j)) * pow(1 - global_m, L - H(i, j));
			}

			D(i, j) = q_ij;
		}

		// Diagonal element
		q_cum = 0;
		for (int64_t j = L - 1; j >= 0; --j)
		{
			q_cum += MutationNeighbours[i][j] * pow(base_u, j + 1) * pow(1 - global_m, L - j - 1);
		}

		D(i, i) = -q_cum;
	}

	QT = MatrixED::Identity(DIM, DIM);
	QT.noalias() += D;
	QT.transposeInPlace();

	// via LU:
	/*
	   Solver fQTinvLU(QT);
	   fQTinv = (fQTinvLU.inverse() + fQTinvLU.inverse().transpose());
	   fQTinv *= 0.5;
	 */

	// via formula:
	/*
	   fQTinv.resize(DIM, DIM);
	   fQTinv = MatrixED::Identity(DIM, DIM);

	   MatrixED TEMP;
	   TEMP = MatrixED::Identity(DIM, DIM);

	   D *= (-1);

	   for (int i = 0; i < 1000; ++i) {
	        TEMP *= D;
	        fQTinv.noalias() += TEMP;
	   }
	 */

	// Check assertions
	bool sameProfiles;
	for (uint64_t i = 0; i < DIM; ++i)
	{
		for (uint64_t j = i; j < DIM; ++j)
		{
			// Symmetry
			if ((D(i, j) != D(j, i)) || (QT(i, j) != QT(j, i)))
			{
				console << "Matrix Q is not symmetrical, aborting\n";
				exit(EXIT_FAILURE);
			}

			if (i != j) {
				// Check mutational profiles for diagonal element
				sameProfiles = true;

				for (uint64_t k = 0; k < L; ++k) {
					sameProfiles = (MutationNeighbours[i][k] == MutationNeighbours[j][k]);

					if (!(sameProfiles))
						break;
				}

				if (sameProfiles) {
					if (verbosity_level >= 3)
						console << "Haplotypes " << i << " and " << j << " have the same mutational pattern\n";

					if ((D(i, i) != D(j, j)) || (QT(i, i) != QT(j, j)))
					{
						console << "but differing self-replication rates, aborting\n";
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	P.resize(DIM);
	P.setIdentity();
	std::random_shuffle(P.indices().data(), P.indices().data() + P.indices().size());

	Q_Random.noalias() = P * QT * P.transpose();

	// fill shuffling vector
	shuffled_index.resize(DIM);
	for (uint32_t i = 0; i < DIM; ++i) {
		shuffled_index[i] = i;
	}

	/*
	   if (verbosity_level >= 3) {
	        console << std::scientific << std::setprecision(6);

	        std::cout << "QT:\n" << QT << '\n';
	        std::cout << "fQTinv:\n" << fQTinv << '\n';
	        std::cout << "fQTinv by LU:\n" << fQTinvLU.inverse() << '\n';
	        std::cout << "Differences in inverses:\n" << fQTinvLU.inverse() - fQTinv << '\n';
	        std::cout << "Shuffling vector:\n" << shuffled_index << '\n';
	   }
	 */
}

