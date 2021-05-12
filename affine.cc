#include <stdio.h>   // printf
#include <algorithm> // shuffle
#include <random>    // rand, default_random_engine
#include <limits>    // numeric_limits
#include <chrono>    // utc timestamps
#include <string>    // probably self-evident

#include "utils.icc"

int main(int arc, char** argv)
{
    // Retrieve Affine parameters from environment vars
    const size_t N_ACTORS = getAffineEnvVar<int>("AFFINE_N_ACTORS");
    const size_t N_EPOCHS = getAffineEnvVar<int>("AFFINE_N_EPOCHS");
    const int EPOCHS_PER_SAVE = getAffineEnvVar<int>("AFFINE_N_EPOCHS_PER_SAVE");
    const int N_LORENZ_POINTS = getAffineEnvVar<int>("AFFINE_N_LORENZ_POINTS");
    const double INIT_WEALTH = getAffineEnvVar<double>("AFFINE_INIT_WEALTH");
    const double TRANSACT_SIZE = getAffineEnvVar<double>("AFFINE_TRANSACT_SIZE");
    const double CHI = getAffineEnvVar<double>("AFFINE_CHI");
    const double ZETA = getAffineEnvVar<double>("AFFINE_ZETA");
    const double KAPPA = getAffineEnvVar<double>("AFFINE_KAPPA");
    // Initialize other Affine-related constants
    const double total_wealth = N_ACTORS*INIT_WEALTH;
    const double loan = (KAPPA/(1. - KAPPA))*INIT_WEALTH;
    // Initialize actors
    double bank[N_ACTORS];
    unsigned int actors[N_ACTORS];
    for (unsigned int actor = 0; actor < N_ACTORS; actor++)
    {
        bank[actor] = INIT_WEALTH;
        actors[actor] = actor;
    }
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
    // Run the simulation
    unsigned seed = getTime();
    XORShiftState* xorshift_state = new XORShiftState(seed);
    for (unsigned int epoch = 0; epoch < N_EPOCHS; epoch++)
    {
        // Shuffle partners
        shuffle(
            &actors[0], &actors[N_ACTORS], 
            std::default_random_engine(seed)
        );
        // Run transactions
        for (unsigned int i = 0; i < N_ACTORS-1; i += 2)
        {
            // Get actor and partner
            unsigned int a = actors[i];
            unsigned int p = actors[i+1];
            // Distribute loans
            bank[a] += loan;
            bank[p] += loan;
            // Determine transaction value
            double transact;
            if (bank[a] < bank[p]) { transact = TRANSACT_SIZE*bank[a]; }
            else { transact = TRANSACT_SIZE*bank[p]; }
            // Determine outcome of biased coin toss
            double prob_minus1 = 0.5*(1 - ZETA*(bank[a] - bank[p])/(INIT_WEALTH));
            bool actor_won = (xorshift64(xorshift_state) > prob_minus1);
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
        // Only compute Lorenz curve at checkpoints or on final epoch
        if (epoch % EPOCHS_PER_SAVE != 0 && epoch != N_EPOCHS-1) { continue; }
        // Print CSV headers
        if (epoch == 0) { printf("epoch,cumulative_wealth,cumulative_actors\n"); }
        // Compute Lorenz curve
        for (unsigned int i = 0; i < N_LORENZ_POINTS; i++)
        {
            double wealth_frac = lorenz_fracs[i]*total_wealth;
            double cumulative_wealth = 0.; // frac of total wealth held by actors w/ wealth < wealth frac
            double cumulative_actors = 0.; // frac of total actors w/ wealth < wealth frac
            for (unsigned int actor = 0; actor < N_ACTORS; actor++)
            {
                if (bank[actor] < wealth_frac) 
                {
                    cumulative_wealth += bank[actor];
                    cumulative_actors += 1.0;
                }
            }
            cumulative_wealth /= total_wealth;
            cumulative_actors /= double(N_ACTORS);
            // Print CSV row
            printf("%d,%f,%f\n", epoch, cumulative_wealth, cumulative_actors);
        }
    }
    return 0;
}
