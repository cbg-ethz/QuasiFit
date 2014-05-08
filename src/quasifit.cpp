#include <quasifit.hpp>

/* Constants */
const uint64_t alph_card = 4;
const EXT_DOUBLE global_m = 3E-5; // global HIV mutation rate
const EXT_DOUBLE base_u = global_m/(alph_card-1); // A -> T mutation rate

void initialise_matrices(const seq_cont& seqs);

int main (int argc, char** argv)
{
	console << "This is " << PACKAGE_STRING << '\n';
	console << "Home page: <" << PACKAGE_URL << ">\n";
	
	/* 1) parse command line */
	parse_arguments(argc, argv);
	
	seq_cont sequences;
	seq_inds indices;
	
	/* 2) process input file */
	VectorED p_MLE;
	load_inputFile(sequences, indices, p_MLE);
	MLE_exists = (Data.minCoeff() != 0);
	
	if ((verbosity_level >= 1) && (MLE_exists))
	{
		for (uint32_t i = 0; i < sequences.size(); ++i)
		{
			console << sequences[i] << '\t' << p_MLE(i) << '\t' << Data(i) << '\n';
		}
	}
	L = sequences[0].length();
	
	if (!(Nreads))
		console << "Sampling from prior\n";
	else
		console << "Sum of samples: " << Nreads << '\n';
	
	// Setup filenames
	uint64_t pos = inputFile.find_last_of('.');
	if (pos != std::string::npos)
	{
		output_m_File = inputFile.substr(0, pos) + "-m.csv";
		output_p_File = inputFile.substr(0, pos) + "-p.csv";
		output_r_File = inputFile.substr(0, pos) + "-r.csv";
		output_ln_File = inputFile.substr(0, pos) + "-ln.txt";
	}
	else
	{
		output_m_File = inputFile + "-m.csv";
		output_p_File = inputFile + "-p.csv";
		output_r_File = inputFile + "-r.csv";
		output_ln_File = inputFile + "-ln.txt";
	}
	
	/* 3) initialize matrices for h(p) */
	// Transition matrix Q
	initialise_matrices(sequences);
	
	/* 4) setup threads and MCMC chains */
	if (GAMMA == -1)
	{
		GAMMA = 1.6829/sqrt(DIM);
	}
	
	if (T > C)
	{
		T = C;
		console << "Reduced number of threads to " << T << '\n';
	}
	
	NChains_per_thread = ceil(1.0*C/T);
	C = NChains_per_thread*T;
	
	no_samples = N/s;
	N = no_samples*s;
	
	console << "Total number of trials (N):             " << N << '\n';
	console << "Total number of chains (C):             " << C << '\n';
	console << "Chains per thread (NChains_per_thread): " << NChains_per_thread << '\n';
	console << "Thinning number (s):                    " << s << '\n';
	console << "Number of samples in total:             " << no_samples << '\n';
	
	matrix_f.resize(DIM, no_samples*C);
	matrix_m.resize(DIM, no_samples*C);
	matrix_s.resize(DIM, no_samples*C);
	matrix_p.resize(DIM, no_samples*C);
	matrix_r.resize(DIM, no_samples*C);
	matrix_ln.resize(3, no_samples*C);
	
	population_f.resize(DIM, C);
	population_m.resize(DIM, C);
	population_s.resize(DIM, C);
	population_p.resize(DIM, C);
	population_r.resize(DIM, C);
	population_r_new.resize(DIM, C);
	population_ln.resize(3, C);
	population_ln_new.resize(3, C);
	
	boost::thread_group* all_threads = new boost::thread_group;
	boost::barrier syn_barrier(T);
	
	counters.resize(C);
	counters.setZero();
	
	acc.resize(C);
	acc.setZero();
	
	// Initialise population
	console << "Initialising population\n";
	inputFile_initial = inputFile + "-initial";
	
	if (load_initial_r_from_file)
	{
		// load initial population from file
		console << "Loading initial population from file\n";
		
		r_loader_from_file();
		
		for (uint64_t t = 0; t < T; ++t)
		{
			all_threads->add_thread(new boost::thread(population_initialiser, t));
		}
		
		if (verbosity_level > 2)
		{
			console << "Log Posterior of initial population: " << population_ln << '\n';
		}
		
		all_threads->join_all();
		delete all_threads;
	}
	else
	{
		// do NOT load initial population from file
		console << "Starting from random position\n";
		
		leave_id = DIM-1;
		
		VectorED r_temp(DIM);
		VectorED p_temp(DIM);
		VectorED m_temp(DIM);
		VectorED f_temp(DIM);
		VectorED s_temp(DIM);
		EXT_DOUBLE temp_logPosterior, temp_logDet, temp_logMult;
		bool reject_outright;
		
		p_temp.noalias() = p_MLE;
		if (!(MLE_exists))
		{
			// use 'pseudocounts' for initialisation
			p_temp += VectorED::Constant(DIM, 1.0/DIM);
			//p_temp /= p_temp.sum();
			
			console << "MLE does not exist, starting from 1/n\n";
		}
		
		// 1) Convert p_temp to r_temp
		convert_from_P_to_R(p_temp, r_temp);
		
		// 2) convert P to M and F
		temp_logPosterior = fitness_space(p_temp, p_temp, f_temp, m_temp, reject_outright, temp_logDet, temp_logMult);
		
		// 4) convert M to S
		convert_from_M_to_S(m_temp, s_temp);
		
		for (uint64_t i = 0; i < C; ++i)
		{
			population_ln(0, i) = temp_logPosterior;
			population_ln(1, i) = temp_logDet;
			population_ln(2, i) = temp_logMult;
			
			population_r.col(i) = r_temp;
			
			population_p.col(i) = p_temp;
			population_f.col(i) = f_temp;
			population_s.col(i) = s_temp;
			population_m.col(i) = m_temp;
		}
	}
	
	console << "Log Posterior of initial population: " << population_ln.row(0) << '\n';
	
	// Initialise timers
	boost::chrono::process_real_cpu_clock::time_point realEnd, realStart = boost::chrono::process_real_cpu_clock::now();
	boost::chrono::process_user_cpu_clock::time_point userEnd, userStart = boost::chrono::process_user_cpu_clock::now();
	boost::chrono::process_system_cpu_clock::time_point systemEnd, systemStart = boost::chrono::process_system_cpu_clock::now();
	
	boost::chrono::high_resolution_clock::time_point timeEnd, timeStart = boost::chrono::high_resolution_clock::now();
	
	// Start MCMC
	if (T > 5)
		console << "Starting MCMC threads\n";
	
	all_threads = new boost::thread_group;
	
	for (uint64_t t = 0; t < T; ++t)
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
	
	secs_t
	lReal(boost::chrono::duration_cast<millisecs_t>(realEnd-realStart)),
	lUser(boost::chrono::duration_cast<millisecs_t>(userEnd-userStart)),
	lSystem(boost::chrono::duration_cast<millisecs_t>(systemEnd-systemStart)),
	lTime(boost::chrono::duration_cast<millisecs_t>(timeEnd-timeStart));
	
	console << std::fixed << std::setprecision(3);
	console << "MCMC ran for:\n";
	console << "     Total Time: " << lTime.count() << "s\n";
	console << "  Real CPU Time: " << lReal.count() << "s\n";
	console << "  User CPU Time: " << lUser.count() << "s\n";
	console << "System CPU Time: " << lSystem.count() << "s\n";
	
	console << "Statistics:\n";
	//console << "Accepted:        " << acc << '\n';
	console << "Acceptance rate: " << static_cast<NORMAL_DOUBLE>(acc.sum())/(N*C)*100 << "%\n";
	
	/* Output the samples */
	std::ofstream output;
	
	if (skip_writing)
	{
		console << "Skipped writing to file\n";
	}
	else
	{
		if (no_headers)
		{
			console << "Writing no headers to output files\n";
		}
		
		Eigen::IOFormat StdFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n", "", "", "", "");
		
		/* fitness samples */
		/*
		 output.open(output_f_File.c_str());
		 output << std::fixed << std::setprecision(20);
		 
		 // header:
		 if (!(no_headers))
		 {
		 output << sequences[indices[0]];
		 for (uint64_t j = 1; j < indices.size(); ++j)
		 {
		 output << delim << sequences[indices[j]];
		 }
		 output << '\n';
		 }
		 
		 // actual data:
		 for (uint64_t i = 0; i < no_samples*C; ++i)
		 {
		 output << matrix_f(indices[0],i);
		 for (uint64_t j = 1; j < indices.size(); ++j)
		 {
		 output << delim << matrix_f(indices[j],i);
		 }
		 output << '\n';
		 }
		 
		 output << matrix_f.transpose().format(StdFormat);
		 output.close();
		 */
		
		/* fitness manifold samples */
		output.open(output_m_File.c_str());
		output << std::fixed << std::setprecision(20);
		
		// header:
		if (!(no_headers))
		{
			output << sequences[indices[0]];
			for (uint64_t j = 1; j < indices.size(); ++j)
			{
				output << delim << sequences[indices[j]];
			}
			output << '\n';
		}
		
		// actual data:
		for (uint64_t i = 0; i < no_samples*C; ++i)
		{
			output << matrix_m(indices[0],i);
			for (uint64_t j = 1; j < indices.size(); ++j)
			{
				output << delim << matrix_m(indices[j],i);
			}
			output << '\n';
		}
		output.close();
		
		/* selection samples */
		/*
		 output.open(output_s_File.c_str());
		 output << std::scientific << std::setprecision(8);
		 
		 // header:
		 if (!(no_headers))
		 {
		 output << sequences[indices[0]];
		 for (uint64_t j = 1; j < indices.size(); ++j)
		 {
		 output << delim << sequences[indices[j]];
		 }
		 output << '\n';
		 }
		 
		 // actual data:
		 for (uint64_t i = 0; i < no_samples*C; ++i)
		 {
		 output << matrix_s(indices[0],i);
		 for (uint64_t j = 1; j < indices.size(); ++j)
		 {
		 output << delim << matrix_s(indices[j],i);
		 }
		 output << '\n';
		 }
		 output.close();
		 */
		
		/* probability samples */
		output.open(output_p_File.c_str());
		output << std::scientific << std::setprecision(10);
		
		// header:
		if (!(no_headers))
		{
			output << sequences[indices[0]];
			for (uint64_t j = 1; j < indices.size(); ++j)
			{
				output << delim << sequences[indices[j]];
			}
			output << '\n';
		}
		
		// actual data:
		for (uint64_t i = 0; i < no_samples*C; ++i)
		{
			output << matrix_p(indices[0],i);
			for (uint64_t j = 1; j < indices.size(); ++j)
			{
				output << delim << matrix_p(indices[j],i);
			}
			output << '\n';
		}
		output.close();
		
		/* euclidean + logPosterior samples */
		output.open(output_r_File.c_str());
		output << std::scientific << std::setprecision(8);
		
		std::ofstream outputln;
		outputln.open(output_ln_File.c_str());
		outputln << std::scientific << std::setprecision(8);
		
		// header:
		outputln << "LogPost,LogDet,LogMult\n";
		// actual data:
		for (uint64_t i = 0; i < no_samples*C; ++i)
		{
			output << matrix_r(0,i);
			for (uint64_t j = 1; j < DIM; ++j)
			{
				output << delim << matrix_r(j, i);
			}
			output << '\n';
			
			outputln << matrix_ln(0, i) << ',' << matrix_ln(1, i) << ',' << matrix_ln(2, i) << '\n';
		}
		output.close();
		outputln.close();
	}
}

void initialise_matrices(const seq_cont& seqs) {
	MatrixID H(DIM, DIM);
	std::vector<VectorID> MutationNeighbours(DIM);
	VectorID temp_Profile(L);
	
	// create Hamming distance matrix
	for (uint64_t i = 0; i < DIM; ++i)
	{
		temp_Profile.setZero();
		
		for (uint64_t j = 0; j < DIM; ++j)
		{
			H(i, j) = hamming_distance(seqs[i], seqs[j]);
			if (i != j)
				++temp_Profile(H(i, j)-1);
		}
		
		MutationNeighbours[i] = temp_Profile;
	}
	
	if (verbosity_level >= 2)
		console << H << '\n';
	console << "Number of haplotypes: " << DIM << '\n';
	
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
				q_ij = pow(base_u, H(i, j)) * pow(1-global_m, L-H(i, j));
			}
			
			D(i, j) = q_ij;
		}
		
		// Diagonal element
		q_cum = 0;
		for (int64_t j = L-1; j >= 0; --j)
		{
			q_cum += MutationNeighbours[i][j] * pow(base_u, j+1)*pow(1-global_m, L-j-1);
		}
		
		D(i, i) = -q_cum;
	}
	
	QT = MatrixED::Identity(DIM, DIM);
	QT.noalias() += D;
	QT.transposeInPlace();
	
	// via LU:
	Solver fQTinvLU(QT);
	fQTinv = (fQTinvLU.inverse() + fQTinvLU.inverse().transpose());
	fQTinv *= 0.5;
	
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
			eigen_assert(D(i, j) == D(j, i));
			eigen_assert(QT(i, j) == QT(j, i));
			
			if (i != j) {
				// Check mutational profiles for diagonal element
				sameProfiles = true;
				
				for (uint64_t k = 0; k < L; ++k) {
					sameProfiles = (MutationNeighbours[i][k] == MutationNeighbours[j][k]);
					
					if (!(sameProfiles))
						break;
				}
				
				if (sameProfiles) {
					std::cout << "Haplotypes " << i << " and " << j << " have the same mutational pattern\n";
					eigen_assert(D(i, i) == D(j, j));
					eigen_assert(QT(i, i) == QT(j, j));
				}
			}
		}
	}
	
	P.resize(DIM);
	P.setIdentity();
	std::random_shuffle(P.indices().data(), P.indices().data()+P.indices().size());
	
	Q_Random.noalias() = P*QT*P.transpose();
	
	// fill shuffling vector
	shuffled_index.resize(DIM);
	for (uint64_t i = 0; i < DIM; ++i) {
		shuffled_index[i] = i;
	}
	
	if (verbosity_level >= 3) {
		console << std::scientific << std::setprecision(6);
		
		std::cout << "QT:\n" << QT << '\n';
		std::cout << "fQTinv:\n" << fQTinv << '\n';
		std::cout << "fQTinv by LU:\n" << fQTinvLU.inverse() << '\n';
		std::cout << "Differences in inverses:\n" << fQTinvLU.inverse()-fQTinv << '\n';
		std::cout << "Shuffling vector:\n" << shuffled_index << '\n';
	}
}
