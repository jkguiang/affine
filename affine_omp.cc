#include <stdio.h>   // printf
#include <omp.h>     // OpenMP pragmas
#include <algorithm> // shuffle
#include <random>    // rand, default_random_engine
#include <limits>    // numeric_limits
#include <chrono>    // utc timestamps
#include <string>    // probably self-evident

#include "utils.icc"

int main(int arc, char** argv)
{
    // Retrieve Affine parameters from environment vars
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
    double bank[N_ACTORS];
    unsigned int actors[N_ACTORS];
    for (unsigned int actor = 0; actor < N_ACTORS; actor++)
    {
        bank[actor] = INIT_WEALTH;
        actors[actor] = actor;
    }
    double transact;
    double prob_minus1;
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
    double cumulative_wealth[N_LORENZ_POINTS] = {0.};
    double cumulative_actors[N_LORENZ_POINTS] = {0.};
    // Initialize iterating indices
    unsigned int i;
    unsigned int a; // actor
    unsigned int p; // partner
    // Run the simulation
    unsigned seed = Utils::getTime();
    Utils::XORShiftState* xorshift_state;
    for (unsigned int epoch = 0; epoch < N_EPOCHS; epoch++)
    {
        // Shuffle partners
        shuffle(
            &actors[0], &actors[N_ACTORS], 
            std::default_random_engine(seed*epoch)
        );
        #pragma omp parallel \
        shared(INIT_WEALTH, TRANSACT_SIZE, CHI, ZETA, bank, loan, actors) \
        private(i, a, p, transact, prob_minus1, actor_won, xorshift_state)
        {
            xorshift_state = new Utils::XORShiftState(seed*omp_get_thread_num());
            // Run transactions
            #pragma omp for
            for (i = 0; i < N_ACTORS-1; i += 2)
            {
                // Get actor and partner
                a = actors[i];
                p = actors[i+1];
                // Distribute loans
                bank[a] += loan;
                bank[p] += loan;
                // Determine transaction value
                if (bank[a] < bank[p]) { transact = TRANSACT_SIZE*bank[a]; }
                else { transact = TRANSACT_SIZE*bank[p]; }
                // Determine outcome of biased coin toss
                prob_minus1 = 0.5*(1 - ZETA*(bank[a] - bank[p])/(INIT_WEALTH));
                actor_won = (Utils::xorshift64(xorshift_state) > prob_minus1);
                // Apply transaction to actor and partner wealth
                bank[a] += transact*pow(-1, !actor_won);
                bank[p] += transact*pow(-1, actor_won);
                // Redistribute wealth
                bank[a] += (INIT_WEALTH+loan - bank[a])*CHI;
                bank[p] += (INIT_WEALTH+loan - bank[p])*CHI;
                // Collect loans
                bank[a] -= loan;
                bank[p] -= loan;
            }
        }
        // Only compute Lorenz curve at checkpoints or on final epoch
        if (epoch % EPOCHS_PER_SAVE != 0 && epoch != N_EPOCHS-1) { continue; }
        // Print CSV headers
        if (epoch == 0) { printf("epoch,cumulative_wealth,cumulative_actors\n"); }
        // Compute Lorenz curve
        #pragma omp parallel \
        shared(N_ACTORS, total_wealth, bank, lorenz_fracs, cumulative_wealth, cumulative_actors) \
        private(i, a)
        {
            #pragma omp for collapse(2)
            for (i = 0; i < N_LORENZ_POINTS; i++)
            {
                for (a = 0; a < N_ACTORS; a++)
                {
                    if (bank[a] < lorenz_fracs[i]*total_wealth) 
                    {
                        cumulative_wealth[i] += bank[a];
                        cumulative_actors[i] += 1.0;
                    }
                }
            }
        }
        for (i = 0; i < N_LORENZ_POINTS; i++)
        {
            cumulative_wealth[i] /= total_wealth;
            cumulative_actors[i] /= double(N_ACTORS);
            // Print CSV row
            printf("%d,%f,%f\n", epoch, cumulative_wealth[i], cumulative_actors[i]);
            // Reset
            cumulative_wealth[i] = 0.;
            cumulative_actors[i] = 0.;
        }
    }
    return 0;
}
