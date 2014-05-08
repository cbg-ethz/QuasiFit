#include <quasifit.hpp>

boost::mutex mutex;

// MCMC options
uint64_t NChains_per_thread;
uint64_t no_samples;

VectorID counters;
VectorID acc;

VectorID shuffled_index;

uint32_t leave_id;

// Data types
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;

MatrixED QT;
MatrixED Q_Random;

MatrixED sQTinv;
MatrixED fQTinv;
Solver Q_LU;

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
MatrixED population_ln;
MatrixED population_ln_new;

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

EXT_DOUBLE calculate_log_multinomial(const VectorED& p)
{
	if (Nreads)
	{
		// we have data, not sampling from the prior
		EXT_DOUBLE sum = 0;
		for (uint64_t i = 0; i < DIM; ++i)
			sum += Data(shuffled_index[i])*log(p(shuffled_index[i]));
		return sum;
	}
	else
	{
		// NO data, sampling from the prior
		return 0;
	}
}

EXT_DOUBLE calculate_log_determinant(const VectorED& f)
{
	EXT_DOUBLE sum = 0;
	
	Solver LU;
	MatrixED tempM = MatrixED::Identity(DIM,DIM);
	MatrixED diagF = f.asDiagonal();
	tempM -= Q_Random * (P * diagF * P.transpose());
	LU.compute(tempM);
	
	/*
	VectorED TempV2 = fQTinv * p;
	MatrixED tempM2 = p.asDiagonal().inverse();
	tempM2 *= TempV2.asDiagonal();
	tempM -= tempM2;
	*/
	
	for (uint64_t i = 0; i < DIM; ++i)
		sum += log(fabs(LU.matrixLU()(shuffled_index[i], shuffled_index[i])));
	
	return sum;
}

EXT_DOUBLE convert_from_probability_to_fitness__manifold(const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, VectorED& M_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult)
{
	logMult = 0;
	logDet = 0;
	
	Solver LU;
	VectorED p_prime = P*P_vector_unnorm;
	VectorED z = Q_Random.partialPivLu().solve(p_prime);
	M_vector = P_vector_unnorm.asDiagonal().inverse() * (P.transpose() * z);
	
	//M_vector.noalias() = P_vector_unnorm.asDiagonal().inverse() * fQTinv * P_vector_unnorm;
	
	if (M_vector.minCoeff() <= 0)
	{
		reject = true;
	}
	else
	{
		reject = false;
		logMult = calculate_log_multinomial(P_vector);
		logDet = calculate_log_determinant(M_vector);
	}
	
	return logMult + logDet;
}

EXT_DOUBLE (*fitness_space) (const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, VectorED& M_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult) = convert_from_probability_to_fitness__manifold;

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
		
		//sum += P_vector_unnorm(i);
	}
	
	VectorED p_t = P_vector_unnorm;
	std::sort(p_t.data(), p_t.data() + p_t.size());
	//std::cout << "p_t:" << '\t' << p_t << '\n';
	for (uint32_t i = 0; i < DIM; ++i)
	{
		sum += p_t(i);
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
void probability_sampler(boost::barrier& syn_barrier, uint64_t thread_no)
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
		console << "Probability sampler\n";
		console << "gamma: " << GAMMA << '\n';
		console << "b:     " << B << '\n';
		if (!(randomise_leave_out))
			console << "NOT randomising position of leave-one-out in r vector\n";
	}
	
	EXT_DOUBLE temp_logPosterior, temp_logDet, temp_logMult;
	
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
			temp_logPosterior = fitness_space(p_new, p_new_unnorm, f_new, m_new, reject_outright, temp_logDet, temp_logMult);
			
			if (verbosity_level >= 4)
			{
				mutex.lock();
				
				std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
				std::cout << "r_new: " << r_new.transpose() << '\n';
				std::cout << "m_new: " << m_new.transpose() << '\n';
				std::cout << "p_new: " << p_new.transpose() << '\n';
				
				mutex.unlock();
			}
			
			if ((!(reject_outright)) && ((temp_logPosterior - population_ln(0, j)) > log(gsl_ran_flat(r, 0, 1))))
			{
				// accept proposal
				// 4) convert M to S
				convert_from_M_to_S(m_new, s_new);
				
				population_ln_new(0, j) = temp_logPosterior;
				population_ln_new(1, j) = temp_logDet;
				population_ln_new(2, j) = temp_logMult;
				population_r_new.col(j) = r_new;
				
				population_p.col(j) = p_new;
				population_f.col(j) = f_new;
				population_s.col(j) = s_new;
				population_m.col(j) = m_new;
				
				++acc(j);
				
				if ((verbosity_level >= 4) && (m_new(1) < 1))
				{
					mutex.lock();
					std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
					std::cout << "r_new: " << r_new.transpose() << '\n';
					std::cout << "m_new: " << m_new.transpose() << '\n';
					std::cout << "p_new: " << p_new.transpose() << '\n';
					mutex.unlock();
				}
			}
			else
			{
				// reject proposal
				population_ln_new(0, j) = population_ln(0, j);
				population_ln_new(1, j) = population_ln(1, j);
				population_ln_new(2, j) = population_ln(2, j);
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
				matrix_ln.block(0, i/s*C, 3, C) = population_ln;
			}
			
			if (randomise_leave_out) {
				leave_id = gsl_rng_uniform_int(r, DIM);
				gsl_ran_shuffle(r, shuffled_index.data(), DIM, sizeof(uint64_t));
				
				// shuffle permutation matrix
				std::random_shuffle(P.indices().data(), P.indices().data()+P.indices().size());
				Q_Random.noalias() = P*QT*P.transpose();
				//std::cout << "Shuffling vector:\n" << shuffled_index << '\n';
			}
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

void (*MCMC)(boost::barrier& syn_barrier, uint64_t thread_no) = probability_sampler;

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
	EXT_DOUBLE temp_logPosterior, temp_logDet, temp_logMult;
	bool reject_outright;
	
	for (int j = C*thread_no/T; j < C*(thread_no+1)/T; ++j)
	{
		// 1) convert R to P + unnormalised P
		convert_from_R_to_P(population_r.col(j), p_new_unnorm, p_new);
		
		// 2) convert P to M and F
		temp_logPosterior = fitness_space(p_new, p_new_unnorm, f_new, m_new, reject_outright, temp_logDet, temp_logMult);
		
		// 3) convert M to S
		convert_from_M_to_S(m_new, s_new);
		
		population_ln(0, j) = temp_logPosterior;
		population_ln(1, j) = temp_logDet;
		population_ln(2, j) = temp_logMult;
		
		population_p.col(j) = p_new;
		population_f.col(j) = f_new;
		population_s.col(j) = s_new;
		population_m.col(j) = m_new;
	}
	
	gsl_rng_free(r);
}
