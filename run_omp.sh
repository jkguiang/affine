if [[ "${1}" != "" ]]; then
    mkdir -p outputs
    # Print a friendly message
    echo "Running with ${1} threads"
    # Set the number of threads
    export OMP_NUM_THREADS=${1}
    # Set Affine environment variables
    export AFFINE_N_EPOCHS=2000
    export AFFINE_N_ACTORS=5000
    export AFFINE_N_EPOCHS_PER_SAVE=10
    export AFFINE_N_LORENZ_POINTS=100
    export AFFINE_INIT_WEALTH=100.0
    export AFFINE_TRANSACT_SIZE=0.25
    export AFFINE_CHI=0.036
    export AFFINE_KAPPA=0.058
    export AFFINE_ZETA=0.05
    # Set module environment
    module purge
    module load slurm
    module load cpu
    module load intel
    module load intel-mpi
    # Compile code
    icpc -qopenmp -std=c++11 -o affine_omp affine_omp.cc
    # Run code
    ./affine_omp &> outputs/omp_lorenz_curves.csv
    # Compress output
    gzip outputs/omp_lorenz_curves.csv

    echo "Done."
else
    echo "Usage: ./run.sh [N_THREADS]"
fi
