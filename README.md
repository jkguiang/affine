# Running on Expanse
This code was written to run on the SDSC Expanse system, so instructions to do so are included below:
1. Clone this repository
```
git clone https://github.com/jkguiang/affine.git
```
2. Obtain an interactive CPU node
```
srun --partition=debug --qos=debug-normal --pty --account=<your account> \
     --nodes=1 --ntasks-per-node=128 --mem=248 -t 00:30:00 --wait=0 --export=ALL /bin/bash
```
3. Run the ```run.sh``` script
    - Running the serial code: `./run.sh --n_reps=10 --sim=usa`
    - Running the OpenMP code: `./run.sh --openmp --n_threads=4 --n_reps=10 --sim=usa`
    - Note: there is also a "debug" simulation that is much smaller and easier to debug.
