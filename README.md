# Individual-based simulations of Fisherian sexual selection
An individual-based model of Fisherian sexual selection. See [Pomiankowski et al (1991) Evolution]( https://doi.org/10.2307/2409889) allowing for costly preferences and mutation biases 

## How to compile
Open a Unix shell (either on a Linux system, MSYS2 on Windows or the terminal on Mac OS X)

    cd fisher_sexsel/src/ibm
    make

## How to run
    cd fisher_sexsel/src/ibm
    ./run_single_simulation.sh

## Vary parameters
Edit the file `run_single_simulation.sh` and look at the function `initArguments()` in the c++ file.
to check which parameter is which.

