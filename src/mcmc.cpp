#include <quasifit.hpp>

// MCMC options
boost::mutex mutex;

uint32_t leave_id;
uint32_t NChains_per_thread;
uint64_t no_samples;

VectorID counters;
VectorID acc;
VectorID shuffled_index;

// Data types
Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P;

MatrixED QT;
MatrixED Q_Random;

// MatrixED sQTinv;
// MatrixED fQTinv;
// Solver Q_LU;

VectorED p_MLE;

MatrixED matrix_f;
MatrixED matrix_p;
MatrixED matrix_r;
MatrixED matrix_ln;

MatrixED population_f;
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
		console << "Could not open /dev/urandom for generating a seed\n";
		exit(EXIT_FAILURE);
	}
	return random_seed;
}

void generate_sample(gsl_rng* r, VectorED& new_sample, const uint32_t chain_id)
{
	uint32_t R1, R2;

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

	new_sample = population_r.col(chain_id) + GAMMA * (population_r.col(R1) - population_r.col(R2));

	// Add Gaussian noise
	for (uint32_t i = 0; i < DIM; ++i)
	{
		if (i != leave_id)
			new_sample(i) += gsl_ran_gaussian(r, B);
	}
}

EXT_DOUBLE calculate_log_multinomial(const VectorED& p)
{
	if (Nreads)
	{
		// we have data, not sampling from the prior
		EXT_DOUBLE sum = 0;
		for (uint32_t i = 0; i < DIM; ++i)
			sum += Data(shuffled_index[i]) * log(p(shuffled_index[i]));
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

	for (uint32_t i = 0; i < DIM; ++i)
		sum += log(fabs(LU.matrixLU() (shuffled_index[i], shuffled_index[i])));

	return sum;
}

void convert_from_P_to_F(const VectorED& p, VectorED& f)
{
	// Solver LU;
	VectorED p_prime = P * p;
	VectorED z = Q_Random.partialPivLu().solve(p_prime);
	f = p.asDiagonal().inverse() * (P.transpose() * z);
}

EXT_DOUBLE convert_from_probability_to_fitness__manifold(const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult)
{
	logMult = 0;
	logDet = 0;

	convert_from_P_to_F(P_vector_unnorm, F_vector);

	if (F_vector.minCoeff() <= 0)
	{
		reject = true;
	}
	else
	{
		reject = false;
		logMult = calculate_log_multinomial(P_vector);
		logDet = calculate_log_determinant(F_vector);
	}

	return logMult + logDet;
}

EXT_DOUBLE (* fitness_space) (const VectorED& P_vector, const VectorED& P_vector_unnorm, VectorED& F_vector, bool& reject, EXT_DOUBLE& logDet, EXT_DOUBLE& logMult) = convert_from_probability_to_fitness__manifold;

void convert_from_R_to_P(const VectorED& R_vector, VectorED& P_vector_unnorm, VectorED& P_vector)
{
	EXT_DOUBLE sum = 0;

	for (uint32_t i = 0; i < DIM; ++i)
	{
		if (i == leave_id)
			P_vector_unnorm(i) = 1;
		else
			P_vector_unnorm(i) = exp(R_vector(i));
	}

	VectorED p_t = P_vector_unnorm;
	std::sort(p_t.data(), p_t.data() + p_t.size());

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

// the actual heart of the sampler
void probability_sampler(boost::barrier& syn_barrier, uint32_t thread_no)
{
	gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
	// gsl_rng* r = gsl_rng_alloc(gsl_rng_ranlxd2);
	// gsl_rng* r = gsl_rng_alloc(gsl_rng_taus);
	uint64_t seed;
	seed = random_seed();
	gsl_rng_set(r, seed);

	bool reject_outright;
	EXT_DOUBLE temp_logPosterior, temp_logDet, temp_logMult;
	VectorED f_new(DIM), p_new(DIM), p_new_unnorm(DIM), r_new(DIM);

	for (uint64_t i = 0; i < N; ++i)
	{
		// Perform N MCMC trials ...
		for (uint64_t j = thread_no * NChains_per_thread; j < (thread_no + 1) * NChains_per_thread; ++j)
		{
			// ... for NChains_per_thread Chains
			++counters(j);

			// 1) generate new R
			generate_sample(r, r_new, j);

			// 2) convert R to P + unnormalised P
			convert_from_R_to_P(r_new, p_new_unnorm, p_new);

			// 3) convert P to F
			temp_logPosterior = fitness_space(p_new, p_new_unnorm, f_new, reject_outright, temp_logDet, temp_logMult);

			if (verbosity_level >= 4)
			{
				mutex.lock();
				std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
				std::cout << "f_new: " << f_new.transpose() << '\n';
				std::cout << "p_new: " << p_new.transpose() << '\n';
				std::cout << "r_new: " << r_new.transpose() << '\n';
				mutex.unlock();
			}

			if ((!(reject_outright)) && ((temp_logPosterior - population_ln(0, j)) > log(gsl_ran_flat(r, 0, 1))))
			{
				// accept proposal
				population_ln_new(0, j) = temp_logPosterior;
				population_ln_new(1, j) = temp_logDet;
				population_ln_new(2, j) = temp_logMult;
				population_r_new.col(j) = r_new;

				population_f.col(j) = f_new;
				population_p.col(j) = p_new;

				++acc(j);

				if ((verbosity_level >= 4) && (f_new(1) < 1))
				{
					mutex.lock();
					std::cout << "LogPosterior: " << log(temp_logPosterior) << '\n';
					std::cout << "f_new: " << f_new.transpose() << '\n';
					std::cout << "p_new: " << p_new.transpose() << '\n';
					std::cout << "r_new: " << r_new.transpose() << '\n';
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
				matrix_f.block(0, i / s * C, DIM, C) = population_f;
				matrix_p.block(0, i / s * C, DIM, C) = population_p;
				matrix_r.block(0, i / s * C, DIM, C) = population_r;
				matrix_ln.block(0, i / s * C, 3, C) = population_ln;
			}

			if (randomise_leave_out) {
				// new leave-one-out-index
				leave_id = gsl_rng_uniform_int(r, DIM);

				// new permutation matrix and indices
				gsl_ran_shuffle(r, shuffled_index.data(), DIM, sizeof(uint32_t));
				// std::cout << "Shuffling vector:\n" << shuffled_index << '\n';
				std::random_shuffle(P.indices().data(), P.indices().data() + P.indices().size());
				Q_Random.noalias() = P * QT * P.transpose();
			}
		}
		syn_barrier.wait();

		// regenerate r population
		if (randomise_leave_out)
		{
			for (uint64_t j = thread_no * NChains_per_thread; j < (thread_no + 1) * NChains_per_thread; ++j)
			{
				convert_from_P_to_R(population_p.col(j), r_new);
				population_r.col(j) = r_new;
			}
			syn_barrier.wait();
		}
	}

	gsl_rng_free(r);
}

void (* MCMC)(boost::barrier& syn_barrier, uint32_t thread_no) = probability_sampler;

void population_initialiser()
{
	console << "Initialising population\n";

	EXT_DOUBLE temp_logPosterior, temp_logDet, temp_logMult;
	VectorED f_temp(DIM);
	VectorED r_temp(DIM);
	p_MLE.resize(DIM);

	if (Data.minCoeff() == 0)
	{
		MLE_exists = false;
	}
	else
	{
		p_MLE = Data.cast<EXT_DOUBLE>() / Data.sum();
		convert_from_P_to_F(p_MLE, f_temp);

		if (f_temp.minCoeff() <= 0)
		{
			MLE_exists = false;
		}
		else
		{
			MLE_exists = true;
		}
	}

	leave_id = DIM - 1;

	/* heuristic for finding a good starting point */
	if (!(MLE_exists))
	{
		// we have to find a good starting point for the sampler, close to the MLE were the space closed
		if (Nreads)
		{
			VectorED p_freq(DIM);
			VectorED p_lower(DIM);

			p_freq = Data.cast<EXT_DOUBLE>() / Data.sum();

			for (uint32_t i = 0; i < DIM; ++i)
			{
				p_lower(i) = closest_observed_neighbour(i) * global_m;
			}

			uint32_t i = 0;
			do
			{
				++i;
				p_MLE = p_freq + i * p_lower;
				p_MLE /= p_MLE.sum();

				convert_from_P_to_F(p_MLE, f_temp);
			}
			while (f_temp.minCoeff() <= 0);
		}
		else
		{
			// prior
			p_MLE = VectorED::Constant(DIM, 1.0 / DIM);
			convert_from_P_to_F(p_MLE, f_temp);
		}
	}

	// 1) Convert p_temp to r_temp
	convert_from_P_to_R(p_MLE, r_temp);

	// 2) convert P to F
	temp_logDet = calculate_log_determinant(f_temp);
	temp_logMult = calculate_log_multinomial(p_MLE);
	temp_logPosterior = temp_logDet + temp_logMult;

	for (uint32_t i = 0; i < C; ++i)
	{
		population_ln(0, i) = temp_logPosterior;
		population_ln(1, i) = temp_logDet;
		population_ln(2, i) = temp_logMult;

		population_f.col(i) = f_temp;
		population_p.col(i) = p_MLE;
		population_r.col(i) = r_temp;
	}
}

