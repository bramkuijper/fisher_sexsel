#ifndef _FISHER_HPP_
#define _FISHER_HPP_

// class that sets up the simulation
// on Fisherian sexual selection, initializes
// the population and underlying genetic values
// then runs the simulation and writes the data

#include <random>
#include <fstream>
#include "parameters.hpp"
#include "individual.hpp"

class Fisher
{
    private:

        // a random device needed to initialize
        // the random number sequence (otherwise always
        // the same random number sequence among replicates)
        std::random_device rd;

        // store the seed used to initialize the random
        // number sequence, so that we can 'replay' the 
        // same simulation for debugging, visualisation purposes
        unsigned int seed;

        // the random number generator itself
        std::mt19937 rng_r;

        // object with all the parameters of this simulation
        Parameters par;

        // the file to write statistics to 
        std::ofstream data_file;

        // uniform distribution
        std::uniform_real_distribution<double> uniform;

        // population of males, females and offspring
        std::vector <Individual> males{};
        std::vector <Individual> females{};
        std::vector <Individual> offspring{};

        // current generation
        unsigned int generation{0};
        double fraction_females_survive{0.0};
        double fraction_males_survive{0};

        void survival();
        void reproduce();
        unsigned int choose(Individual const &female);
        void write_data();
        void write_parameters();
        void write_data_headers();

    public:
        // initialize the simulation
        Fisher(Parameters const &parms);

        // run the thing
        void run_simulation();
};


#endif
