#include <random>
#include "parameters.hpp"
#include "individual.hpp"

Individual::Individual(Parameters const &parms) :
    t{parms.init_t,parms.init_t},
    p{parms.init_p,parms.init_p}
{}

Individual::Individual(Individual const &other) :
    t{other.t[0],other.t[1]},
    p{other.p[0],other.p[1]}
{}

Individual::Individual(Individual const &mother,
        Individual const &father,
        std::mt19937 &rng_r,
        Parameters const &params) 
{
    // set up bernoulli distribution to 
    // randomly sample parental alleles
    std::bernoulli_distribution random_segregator{0.5};
    std::uniform_real_distribution unif{0.0,1.0};

    t[0] = mother.t[random_segregator(rng_r)];
    t[1] = father.t[random_segregator(rng_r)];

    p[0] = mother.p[random_segregator(rng_r)];
    p[1] = father.p[random_segregator(rng_r)];

    // cycle through alleles
    for (unsigned int allele_idx{0}; allele_idx < 2; 
            ++allele_idx)
    {
        // mutation event of p allele
        if (unif(rng_r) < params.mu_p) 
        {
            // generate the absolute size of the
            // increment by which the allele
            // will change
            double increment_p{params.max_mut_p * unif(rng_r)};

            // decide whether we need to add or subtract
            p[allele_idx] += unif(rng_r) < 0.5 ? -increment_p : increment_p;

            if (params.only_positive && p[allele_idx] < 0.0)
            {
                p[allele_idx] = 0.0;
            }
        }

        if (unif(rng_r) < params.mu_t)
        {
            double increment_t{params.max_mut_t * unif(rng_r)};

            // decide whether we need to add or subtract
            t[allele_idx] += unif(rng_r) < params.biast ? -increment_t : increment_t;

            if (params.only_positive && t[allele_idx] < 0.0)
            {
                t[allele_idx] = 0.0;
            }
        }
    } // end for

} // end birth constructor

void Individual::operator=(Individual const &other)
{
    for (unsigned int allele_idx{0};
            allele_idx < 2;
            ++allele_idx)
    {
        p[allele_idx] = other.p[allele_idx];
        t[allele_idx] = other.t[allele_idx];
    }
} 
