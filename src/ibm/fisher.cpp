#include <cassert>
#include <random>
#include <iostream>
#include "fisher.hpp"


Fisher::Fisher(Parameters const &params) :
    rd{}, // initialize random device (with which to choose a random seed)
    seed{rd()}, // get the random seed 
    rng_r{seed}, // now use the random seed to initialize a random number generator
    par{params}, // copy over the parameters
    data_file{par.file_name}, // initialize the data file to write output to
    uniform{0.0,1.0}, // initialize the uniform distribution 
    males{par.N/2, Individual(par)}, // initialize all females
    females{par.N/2, Individual(par)} // initialize all males
{}

void Fisher::run_simulation()
{
    write_data_headers();
    // have individuals survive dependent on their 
    // ornaments and preferences
    // then choose mates
    for (generation = 0; 
            generation <= par.max_generation; 
            ++generation)
    {
        // write statistics to file
        if (generation % par.data_print_interval == 0)
        {
            write_data();
        }
        // have individuals survive (or not)
        survival();

        // replace the current generation
        reproduce();

    } // end for unsigned int generation
   
    write_parameters();
        
} // end run_simulation()

void Fisher::survival()
{
    double t, p, psurvive; // auxiliary variable

    // male survival
    for (auto male_iterator{males.begin()};
            male_iterator != males.end();
            )
    {
        // calculate male phenotype
        t = 0.5 * (male_iterator->t[0] + male_iterator->t[1]);

        psurvive = exp(-par.c * t * t);

        if (uniform(rng_r) > psurvive) // calculate probability individual dies
        {
            // replace current male by male at the end of the stack
            std::swap(*male_iterator, males.back());

            // then remove male at the end of the stack
            males.pop_back();
        }
        else
        {
            ++male_iterator;
        }
    } // end for male_iterator

    for (auto female_iterator{females.begin()};
            female_iterator != females.end();
        )
    {
        p = 0.5 * (female_iterator->p[0] + female_iterator->p[1]);

        psurvive = exp(-par.b * p * p);

        if (uniform(rng_r) > psurvive) // calculate probability individual dies
        {
            // move female that is about to be removed to the end
            // of the stack and replace by final elmt, survival
            // of which will be evaluated during the next iteration
            std::swap(*female_iterator, females.back());

            // remove final elmt
            females.pop_back();
        }
        else
        {
            ++female_iterator;
        }
    } // end for female_iterator

    fraction_females_survive = static_cast<double>(females.size()) / (par.N/2);
    fraction_males_survive = static_cast<double>(males.size()) / (par.N/2);
} // end survival()

unsigned int Fisher::choose(Individual const &female)
{
    // express a female's preference
    double p{0.5*(female.p[0] + female.p[1])};

    /// aux variable to store male index
    unsigned int current_male_idx;

    // aux variable to store male ornament
    double t, male_score;

    // store a vector with scores of each male
    std::vector <double> male_scores;
    // store a vector with indices of each male
    std::vector <unsigned int> male_indices;

    // sample distribution of males
    std::uniform_int_distribution<unsigned int> 
        male_sampler{0, static_cast<unsigned int>(males.size() - 1)};

    for (unsigned int male_idx{0}; 
            male_idx < par.n_males_sampled;
            ++male_idx)
    {
        current_male_idx = male_sampler(rng_r);

        t = 0.5 * (males[current_male_idx].t[0] + males[current_male_idx].t[1]);

        male_score = exp(par.a * p * t);

        male_scores.push_back(male_score);
        male_indices.push_back(current_male_idx);
    }

    // now make a distribution of all male scores
    std::discrete_distribution <unsigned int> 
        male_score_distribution(male_scores.begin(), male_scores.end());

    unsigned int sampled_male{male_score_distribution(rng_r)};

    assert(sampled_male < par.n_males_sampled);

    // draw a sample from the distribution of male scores
    // this is the male that will be chosen
    return(male_indices[male_score_distribution(rng_r)]);

} // end choose()

void Fisher::reproduce()
{
    // make sure there are no offspring from
    // the previous generation
    offspring.clear();

    if (females.size() < 1 
            ||
            males.size() < 1)
    {
        write_data();
        write_parameters();
        exit(1);
    }

    // random sampler of females 
    std::uniform_int_distribution<unsigned int> 
        female_sampler{0, static_cast<unsigned int>(females.size() - 1)};

    for (unsigned int offspring_idx{0}; 
            offspring_idx < par.N;
            ++offspring_idx)
    {
        unsigned int mother_id{female_sampler(rng_r)};

        // perform mate choice
        unsigned int father_id{choose(females[mother_id])};

        assert(mother_id < females.size());
        assert(father_id < males.size());

            // use birth constructor
        Individual Kid{Individual(
                females[mother_id],
                males[father_id],
                rng_r,
                par)};

        offspring.push_back(Kid);
    } // end for unsigned offspring

    assert(offspring.size() == par.N);

    unsigned int offspring_counter{0};

    females.clear();
    males.clear();

    for (unsigned int female_idx{0};
            female_idx < par.N/2;
            ++female_idx)
    {
        females.push_back(offspring[offspring_counter]);
        ++offspring_counter;
    }
    
    for (unsigned int male_idx{0};
            male_idx < par.N/2;
            ++male_idx)
    {
        assert(offspring_counter < par.N);
        males.push_back(offspring[offspring_counter]);
        ++offspring_counter;
    }
} // end reproduce()

void Fisher::write_data()
{
    double mean_t{0.0};
    double ss_t{0.0};
    double mean_p{0.0};
    double ss_p{0.0};
    double sum_pt{0.0};

    double t,p;

    for (auto male_iterator{males.begin()};
            male_iterator != males.end();
            ++male_iterator)
    {
        t = 0.5 * (male_iterator->t[0] + male_iterator->t[1]);
        p = 0.5 * (male_iterator->p[0] + male_iterator->p[1]);
        mean_t += t;
        mean_p += p;
        ss_t += t*t;
        ss_p += p*p;
        sum_pt += p*t;
    }
    
    for (auto female_iterator{females.begin()};
            female_iterator != females.end();
            ++female_iterator)
    {
        t = 0.5 * (female_iterator->t[0] + female_iterator->t[1]);
        p = 0.5 * (female_iterator->p[0] + female_iterator->p[1]);
        mean_t += t;
        mean_p += p;
        ss_t += t*t;
        ss_p += p*p;
        sum_pt += p*t;
    }

    unsigned int n = static_cast<unsigned int>(males.size() + females.size());

    mean_t /= n;
    mean_p /= n;

    double var_t = ss_t / n - mean_t * mean_t;
    double var_p = ss_p / n - mean_p * mean_p;

    double cov_pt = sum_pt / n - mean_p * mean_t;

    data_file << generation << ";" 
        << mean_t << ";" 
        << mean_p << ";" 
        << var_t << ";" 
        << var_p << ";"
        << cov_pt << ";"
        << fraction_females_survive << ";"
        << fraction_males_survive << std::endl;
} // end write_data()

void Fisher::write_data_headers()
{
    data_file << "generation;mean_t;mean_p;var_t;var_p;cov_pt;frac_female_survive;frac_male_survive"
        << std::endl;
} // end write_data_headers()

void Fisher::write_parameters()
{
    data_file << std::endl << std::endl
        << "seed;" << seed << ";" << std::endl
        << "N;" << par.N << ";" << std::endl
        << "max_generation;" << par.max_generation << ";" << std::endl
        << "n_males_sampled;" << par.n_males_sampled << ";" << std::endl
        << "mu_t;" << par.mu_t << ";" << std::endl
        << "mu_p;" << par.mu_p << ";" << std::endl
        << "max_mut_p;" << par.max_mut_p << ";" << std::endl
        << "max_mut_t;" << par.max_mut_t << ";" << std::endl
        << "biast;" << par.biast << ";" << std::endl
        << "a;" << par.a << ";" << std::endl
        << "b;" << par.b << ";" << std::endl
        << "c;" << par.c << ";" << std::endl
        << "only_positive;" << par.only_positive << ";" << std::endl
        << "bias_negative;" << par.bias_negative << ";" << std::endl
        << "init_t;" << par.init_t << ";" << std::endl
        << "init_p;" << par.init_p << ";" << std::endl;
}
