# Default values
openmp=false
n_threads=1
sim="debug"
args=""

# Parse arguments
for arg in "$@"; do
    key=$(echo $arg | cut -f1 -d=)
    val=$(echo $arg | cut -f2 -d=)   
    case "$key" in
        --openmp) openmp=true;;
        --n_threads) n_threads=${val};;
        --sim) sim=${val};;     
        *) args+="$arg ";;
    esac
done

if [[ "${sim}" == "debug" ]]; then
    export AFFINE_N_EPOCHS=50
    export AFFINE_N_ACTORS=128
    export AFFINE_N_EPOCHS_PER_SAVE=10
    export AFFINE_N_LORENZ_POINTS=100
    export AFFINE_INIT_WEALTH=1000.0
    export AFFINE_TRANSACT_SIZE=0.25
    export AFFINE_CHI=0.036
    export AFFINE_ZETA=0.05
    export AFFINE_KAPPA=0.058
elif [[ "${sim}" == "usa" ]]; then
    export AFFINE_N_EPOCHS=200
    export AFFINE_N_ACTORS=8192
    export AFFINE_N_EPOCHS_PER_SAVE=128
    export AFFINE_N_LORENZ_POINTS=64
    export AFFINE_INIT_WEALTH=100.0
    export AFFINE_TRANSACT_SIZE=0.25
    export AFFINE_CHI=0.036
    export AFFINE_KAPPA=0.058
    export AFFINE_ZETA=0.05
else
    echo "No pre-programmed simulation called '${sim}'."
    exit 1
fi

# Set up output dir
mkdir -p outputs

# Set module environment
module purge
module load slurm
module load cpu
module load intel

output_csv=""
if [[ "$openmp" == true ]]; then
    # Load OpenMP compiler
    module load intel-mpi
    # Set number of threads
    export OMP_NUM_THREADS=${n_threads}
    if [[ "${n_threads}" == "1" ]]; then
        echo "Running OpenMP code with ${n_threads} thread..."
        output_csv=outputs/${sim}_omp_${n_threads}thread_lorenz_curves.csv
    else
        echo "Running OpenMP code with ${n_threads} threads..."
        output_csv=outputs/${sim}_omp_${n_threads}threads_lorenz_curves.csv
    fi
    # Compile code
    icpc -qopenmp -qopt-report=5 -std=c++11 -o affine_omp affine_omp.cc
    # Run code
    ./affine_omp &> ${output_csv}
else
    echo "Running serial code..."
    output_csv=outputs/${sim}_lorenz_curves.csv
    # Compile code
    icpc -std=c++11 -o affine affine.cc
    # Run code
    ./affine &> ${output_csv}
fi

echo "Compressing output..."
# Compress output
gzip -f ${output_csv}
echo "Done. Saved to ${output_csv}.gz"
