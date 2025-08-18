#ifndef _PARAMETERS_HPP_
#define _PARAMETERS_HPP_

#include <string>

class Parameters
{
    public:
        // population size, which should be even number
        unsigned int N{5000};

        // duration of the simulation
        unsigned int max_generation{20000};

        unsigned int n_males_sampled{10};

        // interval of number of generations 
        // at which data is printed
        unsigned int data_print_interval{10};

        // per-generation probabability that t allele mutates
        double mu_t{0.05};
        // per-generation probabability that p allele mutates
        double mu_p{0.05};

        // ranges of mutation in uniform distribution
        double max_mut_p{0.4};
        double max_mut_t{0.4};

        // mutation bias on t
        // 0.5: no bias, > 0.5 mutations which result in smaller t more likely
        double biast{0.5};

        // efficacy of sexual selection (see Pomiankowski et al 1991 Evolution)
        double a{1.0};

        // survival cost of female choice (see Pomiankowski et al 1991 Evolution)
        double b{0.0025};

        // survival cost of male ornamentation
        double c{0.5};

        // initial values:
        double init_t{0.0};

        double init_p{1.0};

        // bounds on character value
        bool only_positive{false};

        std::string file_name{"sim_fisher.txt"};
};

#endif
