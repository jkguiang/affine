#include <stdio.h>   // printf
#include <omp.h>     // OpenMP pragmas
#include <algorithm> // shuffle
#include <cstring>   // memcpy
#include <random>    // rand, default_random_engine
#include <limits>    // numeric_limits
#include <chrono>    // utc timestamps
#include <string>    // probably self-evident

struct XORShiftState {
    uint64_t a;
    XORShiftState(uint64_t _a) { a = _a; }
};

double xorshift64(XORShiftState* state) {
    uint64_t x = state->a;   /* The state must be seeded with a nonzero value. */
    x ^= x >> 12; // a
    x ^= x << 25; // b
    x ^= x >> 27; // c
    state->a = x;
    auto u64val = x*0x2545F4914F6CDD1D;
    double dval = (double)u64val/(double)std::numeric_limits<uint64_t>::max();
    return dval;
}

template<typename Type>
Type getAffineEnvVar(std::string const& key) 
{
    printf("WARNING: No affine variable of this type");
    return NULL; 
}

template <>
int getAffineEnvVar<int>(std::string const& key)
{
    char* val = getenv(key.c_str());
    if (val == NULL)
    {
        printf("WARNING: No affine variable called %s\n", key.c_str());
        return NULL;
    }
    return std::stoi(val, NULL);
}

template <>
double getAffineEnvVar<double>(std::string const& key)
{
    char* val = getenv(key.c_str());
    if (val == NULL)
    {
        printf("WARNING: No affine variable called %s\n", key.c_str());
        return NULL;
    }
    return std::strtod(val, NULL);
}

size_t getTime()
{
    return std::chrono::system_clock::now().time_since_epoch().count();
}

int main(int arc, char** argv)
{
    // Retrieve Affine parameters from environment vars
    const size_t N_ACTORS = getAffineEnvVar<int>("AFFINE_N_ACTORS");
    const size_t N_EPOCHS = getAffineEnvVar<int>("AFFINE_N_EPOCHS");
    const int EPOCHS_PER_SAVE = getAffineEnvVar<int>("AFFINE_N_EPOCHS_PER_SAVE");
    const size_t N_LORENZ_POINTS = getAffineEnvVar<int>("AFFINE_N_LORENZ_POINTS");
    const double INIT_WEALTH = getAffineEnvVar<double>("AFFINE_INIT_WEALTH");
    const double TRANSACT_SIZE = getAffineEnvVar<double>("AFFINE_TRANSACT_SIZE");
    const double CHI = getAffineEnvVar<double>("AFFINE_CHI");
    const double ZETA = getAffineEnvVar<double>("AFFINE_ZETA");
    const double KAPPA = getAffineEnvVar<double>("AFFINE_KAPPA");
    printf("%d, %f, %f, %f\n", N_ACTORS, CHI, ZETA, KAPPA);
    // Initialize other Affine-related constants
    const double total_wealth = N_ACTORS*INIT_WEALTH;
    const double loan = (KAPPA/(1. - KAPPA))*INIT_WEALTH;
    // Initialize actors
    double actor_wealth[N_ACTORS];
    unsigned int actors[N_ACTORS];
    for (unsigned int actor = 0; actor < N_ACTORS; actor++)
    {
        actor_wealth[actor] = INIT_WEALTH;
        actors[actor] = actor;
    }
    // Initialize Lorenz Curve
    double cumulative_wealth[N_LORENZ_POINTS];
    double cumulative_actors[N_LORENZ_POINTS];
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
        for (unsigned int i = 0; i < N_ACTORS-1; i++)
        {
            // Get actor and partner
            unsigned int actor = actors[i];
            unsigned int partner = actors[i+1];
            // Distribute loans
            actor_wealth[actor] += loan;
            actor_wealth[partner] += loan;
            // Retrieve actor and partner wealth
            double a_wealth = actor_wealth[actor];
            double p_wealth = actor_wealth[partner];
            // Determine transaction value
            double transact;
            if (a_wealth < p_wealth) { transact = TRANSACT_SIZE*a_wealth; }
            else { transact = TRANSACT_SIZE*p_wealth; }
            // Determine outcome of biased coin toss
            double prob_minus1 = 0.5*(1 - ZETA*(a_wealth - p_wealth)/(INIT_WEALTH));
            bool actor_won = (xorshift64(xorshift_state) > prob_minus1);
            // Apply transaction to actor and partner wealth
            actor_wealth[actor] += transact*pow(-1, !actor_won);
            actor_wealth[partner] += transact*pow(-1, actor_won);
            // Redistribute wealth
            actor_wealth[actor] += (INIT_WEALTH - actor_wealth[actor])*CHI;
            actor_wealth[partner] += (INIT_WEALTH - actor_wealth[partner])*CHI;
            // Collect loans
            actor_wealth[actor] -= loan;
            actor_wealth[partner] -= loan;
        }
        // Only compute Lorenz curve at checkpoints or on final epoch
        if (epoch % EPOCHS_PER_SAVE != 0 && epoch != N_EPOCHS - 1) { continue; }
        // Compute Lorenz curve
        for (unsigned int i = 0; i < N_LORENZ_POINTS; i++)
        {
            double wealth_frac = pow(2., i)/pow(2., N_LORENZ_POINTS - 1)*total_wealth;
            cumulative_wealth[i] = 0.;
            cumulative_actors[i] = 0.;
            for (unsigned int actor = 0; actor < N_ACTORS; actor++)
            {
                if (actor_wealth[actor] < wealth_frac) 
                {
                    cumulative_wealth[i] += actor_wealth[actor];
                    cumulative_actors[i] += 1.0;
                }
            }
            cumulative_wealth[i] /= total_wealth;
            cumulative_actors[i] /= (double)N_ACTORS;
        }
    }
    return 0;
}
