# Individual-based simulations of Fisherian sexual selection
An individual-based model of Fisherian sexual selection. See [Pomiankowski et al (1991) Evolution]( https://doi.org/10.2307/2409889), allowing for costly preferences and mutation biases 

Open a Unix shell (either on a Linux system, MSYS2 on Windows or the terminal on Mac OS X). For UExeter users, best is to use MobaXterm and move the c++ files to the rstudio server.

## How to compile
Assuming one has downloaded the files to one's home directory (`~`), use the shell to:

    cd ~/fisher_sexsel/src/ibm/
    make

## How to run
    cd ~/fisher_sexsel/src/ibm/
    ./run_single_simulation.sh

this will produce an output file (typically called `output.csv`) with means, variances and covariances. The separator is a semicolon, so when plotting in `R` take care to use `read_delim` with the correct arguments.

## Vary parameters
Edit the file `run_single_simulation.sh` and look at the function `initArguments()` in the c++ file.
to check which parameter is which.

