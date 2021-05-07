if [[ "${1}" != "" ]]; then
    # Print a friendly message
    echo "Running with ${1} threads"

    # Set the number of threads
    export OMP_NUM_THREADS=${1}
    export AFFINE_N_EPOCHS=300
    export AFFINE_N_ACTORS=100
    export AFFINE_N_EPOCHS_PER_SAVE=100
    export AFFINE_N_LORENZ_POINTS=20
    export AFFINE_INIT_WEALTH=1000.0
    export AFFINE_TRANSACT_SIZE=0.25
    export AFFINE_CHI=0.036
    export AFFINE_ZETA=0.05
    export AFFINE_KAPPA=0.058
    # Set module environment
    module purge
    module load slurm
    module load cpu
    module load intel
    module load intel-mpi

    # Compile code
    icpc -qopenmp -std=c++11 -o affine affine.cc

    ./affine

    echo "Done."
else
    echo "Usage: ./run.sh [N_THREADS]"
fi
