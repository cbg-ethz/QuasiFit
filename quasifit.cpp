#include <config.h>
#include "quasifit.hpp"

int main (int argc, char** argv)
{
    parse_arguments(argc, argv);

    seq_cont sequences;
    seq_inds indices;

    /* Process input file */
    VectorED p_MLE;

    load_inputFile(sequences, indices, p_MLE);

    if (verbosity_level >= 1)
    {
        for (uint32_t i = 0; i < sequences.size(); ++i)
        {
            std::cout << sequences[i] << '\t' << p_MLE(i) << '\t' << Data(i) << '\n';
        }
    }
    L = sequences[0].length();
    std::cout << "Sum of samples: " << Data.sum() << '\n';

    if (Data.minCoeff() != 0)
    {
        MLE_exists = true;
    }
    else
    {
        MLE_exists = false;
    }

    /* Setup filenames */
    uint64_t pos = inputFile.find_last_of('.');
    if (pos != std::string::npos)
    {
        output_f_File = inputFile.substr(0, pos) + "-f.csv";

        output_m_File = inputFile.substr(0, pos) + "-m.csv";
        output_s_File = inputFile.substr(0, pos) + "-s.csv";

        output_p_File = inputFile.substr(0, pos) + "-p.csv";
        output_r_File = inputFile.substr(0, pos) + "-r.csv";

        output_ln_File = inputFile.substr(0, pos) + "-ln.txt";
    }
    else
    {
        output_f_File = inputFile + "-f.csv";

        output_m_File = inputFile + "-m.csv";
        output_s_File = inputFile + "-s.csv";

        output_p_File = inputFile + "-p.csv";
        output_r_File = inputFile + "-r.csv";

        output_ln_File = inputFile + "-ln.txt";
    }

    /* Setup hamming distance matrix H */
    MatrixID H(DIM, DIM);
    for (uint64_t i = 0; i < DIM; ++i)
    {
        for (uint64_t j = 0; j < DIM; ++j)
        {
            H(i, j) = hamming_distance(sequences[i], sequences[j]);
        }
    }

    if (verbosity_level >= 2)
        std::cout << H << '\n';

    std::cout << "Number of haplotypes: " << DIM << '\n';

    if (GAMMA == -1)
    {
        GAMMA = 1.6829/sqrt(DIM);
    }

    MatrixED D(DIM, DIM);
    /* Setup transition matrix Q */
    EXT_DOUBLE q_cum, q_ij;

    for (uint64_t i = 0; i < DIM; ++i)
    {
        q_cum = 0;

        for (uint64_t j = 0; j < DIM; ++j)
        {
            if (i == j)
            {
                q_ij = 0;
            }
            else
            {
                q_ij = pow(base_u, H(i, j)) * pow(1-global_m, L-H(i, j));
            }
            q_cum += q_ij;

            D(i, j) = q_ij;
        }

        D(i, i) = -q_cum;
    }
    //std::cout << Q << "\n";
    QT = MatrixED::Identity(DIM, DIM);// + D;
    QT.noalias() += D;

    /* Precalculate inverse matrices */
    QT.transposeInPlace();

    fQTinv.resize(DIM, DIM);
    fQTinv = MatrixED::Identity(DIM, DIM);

    MatrixED TEMP;
    TEMP = MatrixED::Identity(DIM, DIM);

    D *= (-1);

    for (int i = 0; i < 100; ++i)
    {
        TEMP *= D;
        fQTinv.noalias() += TEMP;
    }

    //fQTinv.resize(DIM, DIM);
    //fQTinv = MatrixED::Identity(DIM, DIM);
    //fQTinv += sQTinv;

    //fQTinv = QT.inverse();
    //sQTinv = fQTinv - MatrixED::Identity(DIM, DIM);

    if (verbosity_level >= 3)
    {
        std::cout << std::fixed << std::setprecision(15);

        std::cout << "QT:\n" << QT.cast<long double>() << '\n';
        std::cout << "fQTinv:\n" << fQTinv.cast<long double>() << '\n';
        //std::cout << "sQTinv" << sQTinv.cast<long double>() << '\n';
    }

    if (T > C)
    {
        T = C;
        std::cout << "Reduced number of threads to " << T << '\n';
    }

    NChains_per_thread = ceil(1.0*C/T);
    C = NChains_per_thread*T;

    no_samples = N/s;
    N = no_samples*s;

    std::cout << "Total number of trials (N):             " << N << '\n';
    std::cout << "Total number of chains (C):             " << C << '\n';
    std::cout << "Chains per thread (NChains_per_thread): " << NChains_per_thread << '\n';
    std::cout << "Thinning number (s):                    " << s << '\n';
    std::cout << "Number of samples in total:             " << no_samples << '\n';

    matrix_f.resize(DIM, no_samples*C);
    matrix_m.resize(DIM, no_samples*C);
    matrix_s.resize(DIM, no_samples*C);
    matrix_p.resize(DIM, no_samples*C);
    matrix_r.resize(DIM, no_samples*C);
    matrix_ln.resize(1, no_samples*C);

    population_f.resize(DIM, C);
    population_m.resize(DIM, C);
    population_s.resize(DIM, C);
    population_p.resize(DIM, C);
    population_r.resize(DIM, C);
    population_r_new.resize(DIM, C);
    population_ln.resize(1, C);
    population_ln_new.resize(1, C);

    boost::thread_group* all_threads = new boost::thread_group;
    boost::barrier syn_barrier(T);

    counters.resize(C);
    counters.setZero();

    acc.resize(C);
    acc.setZero();

    // Initialise population
    std::cout << "Initialising population\n";
    inputFile_initial = inputFile + "-initial";

    std::cout /*<< "Log Posterior of initial population: "*/ << population_ln << '\n';

    if (load_initial_r_from_file)
    {
        // load initial population from file
        r_loader_from_file();

        for (uint64_t t = 0; t < T; ++t)
        {
            all_threads->add_thread(new boost::thread(population_initialiser, t));
        }

        if (verbosity_level > 2)
        {
            std::cout << "Log Posterior of initial population: " << population_ln << '\n';
        }

        all_threads->join_all();
        delete all_threads;
    }
    else
    {
        // do NOT load initial population from file
        leave_id = DIM-1;

        VectorED r_temp(DIM);
        VectorED p_temp(DIM);
        VectorED m_temp(DIM);
        VectorED f_temp(DIM);
        VectorED s_temp(DIM);
        EXT_DOUBLE tempLogPosterior;
        bool reject_outright;

        p_temp.noalias() = p_MLE;
        if (!(MLE_exists))
        {
            // use 'pseudocounts' for initialisation
            p_temp += VectorED::Constant(DIM, 5E-5);
            p_temp /= p_temp.sum();
        }

        // 1) Convert p_temp to r_temp
        convert_from_P_to_R(p_temp, r_temp);

        // 2) convert P to M and F
        tempLogPosterior = fitness_space(p_temp, f_temp, m_temp, reject_outright);

        // 4) convert M to S
        convert_from_M_to_S(m_temp, s_temp);

        for (uint64_t i = 0; i < C; ++i)
        {
            population_ln(0, i) = tempLogPosterior;
            population_r.col(i) = r_temp;

            population_p.col(i) = p_temp;
            population_f.col(i) = f_temp;
            population_s.col(i) = s_temp;
            population_m.col(i) = m_temp;
        }
    }

    // Initialise timers
    boost::chrono::process_real_cpu_clock::time_point realEnd, realStart = boost::chrono::process_real_cpu_clock::now();
    boost::chrono::process_user_cpu_clock::time_point userEnd, userStart = boost::chrono::process_user_cpu_clock::now();
    boost::chrono::process_system_cpu_clock::time_point systemEnd, systemStart = boost::chrono::process_system_cpu_clock::now();

    boost::chrono::high_resolution_clock::time_point timeEnd, timeStart = boost::chrono::high_resolution_clock::now();

    // Start MCMC
    if (T > 5)
        console("Starting MCMC threads");

    all_threads = new boost::thread_group;
    //std::cout << "Number of threads in thread_group: " << all_threads->size() << '\n';

    for (uint64_t t = 0; t < T; ++t)
    {
        if (T <= 5)
            std::cout << "Starting thread #" << t << '\n';

        all_threads->add_thread(new boost::thread(MCMC,
                                boost::ref(syn_barrier),
                                t));
    }

    if (verbosity_level >= 1)
        all_threads->add_thread(new boost::thread(counter_display));

    std::cout << "Processing all threads\n";
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

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "MCMC ran for:\n";
    std::cout << "     Total Time: " << lTime.count() << "s\n";
    std::cout << "  Real CPU Time: " << lReal.count() << "s\n";
    std::cout << "  User CPU Time: " << lUser.count() << "s\n";
    std::cout << "System CPU Time: " << lSystem.count() << "s\n";

    std::cout << "Statistics:\n";
    //std::cout << "Accepted:        " << acc << '\n';
    std::cout << "Acceptance rate: " << static_cast<long double>(acc.sum())/(N*C)*100 << "\%\n";

    /* Output the samples */
    std::ofstream output;

    if (skip_writing)
    {
        std::cout << "Skipped writing to file\n";
    }
    else
    {
        if (no_headers)
        {
            std::cout << "Writing no headers to output files\n";
        }

        /* fitness samples */
        output.open(output_f_File.c_str());
        output << std::fixed << std::setprecision(14);

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
            output << static_cast<NORMAL_DOUBLE>(matrix_f(indices[0],i));
            for (uint64_t j = 1; j < indices.size(); ++j)
            {
                output << delim << static_cast<NORMAL_DOUBLE>(matrix_f(indices[j],i));
            }
            output << '\n';
        }
        output.close();

        /* fitness manifold samples */
        output.open(output_m_File.c_str());
        output << std::fixed << std::setprecision(14);

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
            output << static_cast<NORMAL_DOUBLE>(matrix_m(indices[0],i));
            for (uint64_t j = 1; j < indices.size(); ++j)
            {
                output << delim << static_cast<NORMAL_DOUBLE>(matrix_m(indices[j],i));
            }
            output << '\n';
        }
        output.close();

        /* selection samples */
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
            output << static_cast<NORMAL_DOUBLE>(matrix_s(indices[0],i));
            for (uint64_t j = 1; j < indices.size(); ++j)
            {
                output << delim << static_cast<NORMAL_DOUBLE>(matrix_s(indices[j],i));
            }
            output << '\n';
        }
        output.close();

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

        for (uint64_t i = 0; i < no_samples*C; ++i)
        {
            output << matrix_r(0,i);
            for (uint64_t j = 1; j < DIM; ++j)
            {
                output << delim << matrix_r(j, i);
            }
            output << '\n';

            outputln << matrix_ln(0, i) << '\n';
        }
        output.close();
        outputln.close();
    }
}
