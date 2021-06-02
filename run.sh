# Default values
openmp=false
n_threads=1
n_reps=1
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
        --n_reps) n_reps=${val};;
        *) args+="$arg ";;
    esac
done

if [[ "${sim}" == "debug" ]]; then
    export AFFINE_N_EPOCHS=48
    export AFFINE_N_ACTORS=128
    export AFFINE_N_EPOCHS_PER_SAVE=24
    export AFFINE_N_LORENZ_POINTS=100
    export AFFINE_INIT_WEALTH=1000.0
    export AFFINE_TRANSACT_SIZE=0.25
    export AFFINE_CHI=0.036
    export AFFINE_ZETA=0.05
    export AFFINE_KAPPA=0.058
elif [[ "${sim}" == "usa" ]]; then
    export AFFINE_N_EPOCHS=256
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
output_name=""
executable=""
if [[ "$openmp" == true ]]; then
    executable=affine_omp_simple
    echo "Compiling OpenMP code..."
    # Load OpenMP compiler
    module load intel-mpi
    # Set number of threads
    export OMP_NUM_THREADS=${n_threads}
    if [[ "${n_threads}" == "1" ]]; then
        echo "Using ${n_threads} thread..."
        output_name=${sim}_omp_${n_threads}thread
    else
        echo "Using ${n_threads} threads..."
        output_name=${sim}_omp_${n_threads}threads
    fi
    # Compile code
    icpc -qopenmp -qopt-report=5 -std=c++11 -o ${executable} ${executable}.cc
    # Move optimization report to outputs directory
    mv ${executable}.optrpt outputs/${output_name}.optrpt
else
    executable=affine
    echo "Running serial code..."
    output_name=${sim}_serial
    # Compile code
    icpc -std=c++11 -o ${executable} ${executable}.cc
fi

runtime=0
for rep in $(seq 1 ${n_reps}); do 
    echo "Running..."
    # Run code
    start_utc_ns=$(date +%s%N)
    ./${executable} &> outputs/${output_name}_lorenz_curves_rep${rep}.csv
    end_utc_ns=$(date +%s%N)
    runtime=$((${end_utc_ns} - ${start_utc_ns}))
    echo "Finished. (${runtime} ns)"
    # Write runtime to file
    echo "${runtime}" >> outputs/${output_name}_runtimes.txt
done
# Wrap up
echo "Compressing outputs..."
gzip -f outputs/*.csv
