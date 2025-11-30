#!/usr/bin/env Rscript

# script that makes a shell file (bash) 
# in which each line
# is the command to run a single simulation
# for a particular parameter combination

# this script will generate such commands
# for all parameter combinations specified above

# see the main.cpp file for the command line parameters
# that are accepted by the simulation

# bias on t
biast <- c(seq(0.5,0.9,0.1),0.99)

# efficacy of sexual selection
a <- 1.0

# cost of female preference
b <- c(0.0025)

# cost of male ornament
c <- c(0.5)

# initial values of t, p
init_t <- 0.0
init_p <- 1.0

# per-allele, per-generation 
# mutation probabilities
# of p,t 
mu_p <- 0.05
mu_t <- 0.05

# maximum mutational stepsize in p, t
sdmu_p <- 0.4
sdmu_t <- 0.4

# number of replicate simulations for each unique
# parameter combination
nrep <- 20 

# whether traits should only attain positive values
only_positive <- F

# generate unique output files that are not easily
# overwritten when generating (and running) a subsequent
# batch of simulations.
# Generate a date_time stamp as a character, put
# that in the output file name and then end the file name
date_str <- format(Sys.time(), "%d%m%Y_%H%M%S")

# generate names for the output file based on date and time 
# to ensure some sort of uniqueness, to avoiding that output
# files from different runs overwrite eachother
output_file_prefix = paste0("sim_fisher_",date_str)

counter <- 0

exe = "./fisher_sexsel.exe"

batch_file_contents <- ""

for (rep_i in 1:nrep)
{
    for (a_i in a)
    {
        for (b_i in b)
        {
            for (c_i in c)
            {
                for (biast_i in biast)
                {
                    counter <- counter + 1
                    file_name_i <- paste0(output_file_prefix,"_",counter)

                    echo_str <- paste("echo",counter)

                    command_str <- paste(exe,
                                    a_i,
                                    b_i,
                                    c_i,
                                    biast_i,
                                    mu_p,
                                    mu_t,
                                    sdmu_p,
                                    sdmu_t,
                                    init_p,
                                    init_t,
                                    as.numeric(only_positive),
                                    file_name_i)

                    # append to batch file contents
                    batch_file_contents <- paste0(batch_file_contents
                            ,"\n"
                            ,echo_str
                            ,"\n"
                            ,command_str)
                } # end for biast
            } # end for c_i
        } # end for b_i
    } # end for a_i
} # end for rep_i

writeLines(text=batch_file_contents)
