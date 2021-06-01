#include <stdio.h>   // printf
#include <omp.h>     // OpenMP pragmas
#include <algorithm> // shuffle
#include <random>    // rand, default_random_engine
#include <limits>    // numeric_limits
#include <chrono>    // utc timestamps
#include <cstring>   // memcpy
#include <string>    // probably self-evident

#include "utils.icc"

int main(int arc, char** argv)
{
    // Retrieve Affine parameters from environment vars
    const size_t N_THREADS = Utils::getAffineEnvVar<int>("OMP_NUM_THREADS");
    const size_t N_ACTORS = Utils::getAffineEnvVar<int>("AFFINE_N_ACTORS");
    const size_t N_EPOCHS = Utils::getAffineEnvVar<int>("AFFINE_N_EPOCHS");
    const int EPOCHS_PER_SAVE = Utils::getAffineEnvVar<int>("AFFINE_N_EPOCHS_PER_SAVE");
    const int N_LORENZ_POINTS = Utils::getAffineEnvVar<int>("AFFINE_N_LORENZ_POINTS");
    const double INIT_WEALTH = Utils::getAffineEnvVar<double>("AFFINE_INIT_WEALTH");
    const double TRANSACT_SIZE = Utils::getAffineEnvVar<double>("AFFINE_TRANSACT_SIZE");
    const double CHI = Utils::getAffineEnvVar<double>("AFFINE_CHI");
    const double ZETA = Utils::getAffineEnvVar<double>("AFFINE_ZETA");
    const double KAPPA = Utils::getAffineEnvVar<double>("AFFINE_KAPPA");
    // Initialize other Affine-related constants
    const double total_wealth = N_ACTORS*INIT_WEALTH;
    const double loan = (KAPPA/(1. - KAPPA))*INIT_WEALTH;
    // Initialize Affine variables
    double bank_orig[N_ACTORS];
    double bank_next[N_ACTORS];
    unsigned int actors[N_ACTORS];
    for (unsigned int actor = 0; actor < N_ACTORS; actor++)
    {
        bank_orig[actor] = INIT_WEALTH;
        bank_next[actor] = INIT_WEALTH;
        actors[actor] = actor;
    }
    double transact;
    double prob_minus_1;
    bool actor_won;
    // Initialize Lorenz curve variables
    double lorenz_fracs[N_LORENZ_POINTS];
    for (unsigned int i = 0; i < N_LORENZ_POINTS-5; i++)
    {
        lorenz_fracs[i] = (1./N_ACTORS)*i/(N_LORENZ_POINTS-5);
    }
    lorenz_fracs[N_LORENZ_POINTS-1] = 1.0;
    lorenz_fracs[N_LORENZ_POINTS-2] = 0.8;
    lorenz_fracs[N_LORENZ_POINTS-3] = 0.5;
    lorenz_fracs[N_LORENZ_POINTS-4] = 0.25;
    lorenz_fracs[N_LORENZ_POINTS-5] = 0.1;
    double cumulative_wealth = 0.;
    double cumulative_actors = 0.;
    // Initialize iterating indices
    unsigned int i; // generic index
    unsigned int a; // actor
    unsigned int p; // partner
    unsigned int t; // thread
    // Run the simulation
    int seed = Utils::getTime();
    const uint64_t max_uint64 = std::numeric_limits<uint64_t>::max();
    uint64_t xorstates_orig[N_THREADS];
    uint64_t xorstates_next[N_THREADS];
    for (i = 0; i < N_THREADS; i++)
    {
        xorstates_orig[N_THREADS] = seed*(i+1);
    }
    for (unsigned int epoch = 0; epoch < N_EPOCHS; epoch++)
    {
        // Shuffle partners
        shuffle(
            &actors[0], &actors[N_ACTORS], 
            std::default_random_engine(seed)
        );
        // Compute the next timestep
        #pragma omp parallel \
        default(none) \
        shared(INIT_WEALTH, TRANSACT_SIZE, CHI, ZETA, loan, actors, bank_orig, bank_next, max_uint64, xorstates_orig, xorstates_next) \
        private(i, a, p, t, transact, prob_minus_1, actor_won)
        {
            t = omp_get_thread_num();
            xorstates_next[t] = xorstates_orig[t];
            // Run transactions
            #pragma omp simd lastprivate(bank_next, xorstates_next)
            for (i = 0; i < N_ACTORS-1; i += 2)
            {
                // Get actor and partner
                a = actors[i];
                p = actors[i+1];
                // Distribute loans
                bank_next[a] = bank_orig[a] + loan;
                bank_next[p] = bank_orig[p] + loan;
                // Determine transaction value
                if (bank_next[a] < bank_next[p])
                {
                    transact = TRANSACT_SIZE*bank_next[a];
                }
                else
                {
                    transact = TRANSACT_SIZE*bank_next[p];
                }
                // Determine outcome of biased coin toss
                prob_minus_1 = 0.5*(1-ZETA*(bank_next[a]-bank_next[p])/(INIT_WEALTH));
                // Toss a coin (using XOR-shift RNG)
                xorstates_next[t] ^= xorstates_next[t] >> 12; // a
                xorstates_next[t] ^= xorstates_next[t] << 25; // b
                xorstates_next[t] ^= xorstates_next[t] >> 27; // c
                auto u64val = xorstates_next[t]*0x2545F4914F6CDD1D;
                actor_won = (double(u64val)/double(max_uint64) > prob_minus_1);
                // Apply transaction to actor and partner wealth
                bank_next[a] += transact*pow(-1, !actor_won);
                bank_next[p] += transact*pow(-1, actor_won);
                // Redistribute wealth
                bank_next[a] += (INIT_WEALTH+loan - bank_next[a])*CHI;
                bank_next[p] += (INIT_WEALTH+loan - bank_next[p])*CHI;
                // Collect loans
                bank_next[a] -= loan;
                bank_next[p] -= loan;
            }
        }
        // Set original arrays equal to the next timestep computed above
        std::memcpy(bank_orig, bank_next, sizeof(double)*N_ACTORS);
        std::memcpy(xorstates_orig, xorstates_next, sizeof(uint64_t)*N_THREADS);
        // Only compute Lorenz curve at checkpoints or on final epoch
        if (epoch % EPOCHS_PER_SAVE != 0 && epoch != N_EPOCHS-1) { continue; }
        // Print CSV headers
        if (epoch == 0) { printf("epoch,cumulative_wealth,cumulative_actors\n"); }
        // Compute Lorenz curve
        #pragma omp parallel \
        default(none) \
        shared(N_ACTORS, N_LORENZ_POINTS, total_wealth, bank_next, lorenz_fracs, epoch) \
        private(i, a, cumulative_wealth, cumulative_actors)
        {
            #pragma omp for
            for (i = 0; i < N_LORENZ_POINTS; i++)
            {
                cumulative_wealth = 0.;
                cumulative_actors = 0.;
                for (a = 0; a < N_ACTORS; a++)
                {
                    if (bank_next[a] < lorenz_fracs[i]*total_wealth) 
                    {
                        cumulative_wealth += bank_next[a];
                        cumulative_actors += 1.;
                    }
                }
                cumulative_wealth /= total_wealth;
                cumulative_actors /= double(N_ACTORS);
                // Print CSV row
                printf("%d,%f,%f\n", epoch, cumulative_wealth, cumulative_actors);
            }
        }
    }
    return 0;
}
