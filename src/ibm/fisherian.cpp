//      Fisherian sexual selection for gonochorists
//
//      Bram Kuijper, Lukas Scharer and Ido Pen
//
//      This work is licensed under a Creative Commons 
//      Attribution-NonCommercial-ShareAlike 4.0 
//      International License. 
//      http://creativecommons.org/licenses/by-nc-sa/4.0/

// load some libraries that we need
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <random>



// random number generators
//
// we need something to 'seed' our random
// number generators (i.e., which sequence 
// of random numbers will we pick)
std::random_device rd; 
unsigned seed = rd();

// the actual random number generator
std::mt19937 rng_r(seed);

// allocate a uniform [0,1] random number distribution
std::uniform_real_distribution <double> uniform(0.0,1.0);

// allocate a bernoulli distribution for allele 
// segregration during inheritance
std::bernoulli_distribution segregator(0.5);



// parameters (some of these values are changed upon initialization
// through the command line, see the argv and argc arrays in the 
// init function down below)
const int N = 5000; // population size
const int N_mate_sample = 10; // number of mates sampled
const int clutch_size = 10; // number of offspring produced

double init_t = 0.0; // initial value for ornament
double init_p = 0.0; // initial value for preference
double a = 1.0; // efficacy of sexual selection
double b = 0.5; // cost of preference 
double c = 0.5; // cost of trait
double biast = 0; // mutation bias 

double mu_p 	  = 0.05;            // mutation rate preference
double mu_t 	  = 0.05;            // mutation rate ornament
double sdmu_p         = 0.4;			 // standard deviation mutation stepsize
double sdmu_t         = 0.4;			 // standard deviation mutation stepsize
                                         //
const double NumGen = 150000; // number of generations
const int skip = 10; // n generations interval before data is printed
double sexlimp = 0; // degree of sex-limited expression in p,t
double sexlimt = 0;

// default name of the output file
std::string file_name = "output.csv";

// preference function:
// 0: open-ended preferences
// 1: absolute
// 2: relative
int pref = 0;

// statistics variables
double meanornsurv = 0;

int popsize = N; // population size between 
bool do_stats = 0;

int generation = 0;
int Nfemales = N / 2, Nmales = N / 2;
int msurvivors = 0;
int fsurvivors = 0;

// TODO: change this to account 4 separate-sex organisms
int father_eggs[N];
int mother_eggs[N];


// the components of an actual individual
struct Individual
{
	double t[2]; // diploid, additive loci for t,p
	double p[2];
    double t_expr; // and store their expressed values
    double p_expr;
	// amount to be allocated to male vs. female function
};

// generate the population
typedef Individual Population[N];

// declare various arrays of N individuals 
Population Females, Males, FemaleSurvivors, MaleSurvivors;

// make an array to store the index values of the parents
// to generate offspring (rather than to make an array with lots of
// offspring an array with lots of indices is cheaper).
int Parents[N*clutch_size][2]; 


// declaring parameters and arrays done, now let's move
// on to declaring functions

// function which obtains arguments from the command line
// for definitions of the various parameters see top of the file
void initArguments(int argc, char *argv[])
{
	a = std::stod(argv[1]);
	b = std::stod(argv[2]);
	c = std::stod(argv[3]);
	biast = std::stod(argv[4]);
	mu_p = std::stod(argv[5]);
	mu_t = std::stod(argv[6]);
	sdmu_p = std::stod(argv[7]);
	sdmu_t = std::stod(argv[8]);
	sexlimp = std::stod(argv[9]);
	sexlimt = std::stod(argv[10]);
    pref = std::stoi(argv[11]);
    file_name = argv[12];
} // end initArguments

// mutation function:
// mutate an allele G given a particular mutation rate mu
// a standard deviation of the mutational effect distribution
// of sdmu and a certain mutational bias
void mutate(double &G, double mu, double sdmu, double bias=0.0)
{
    // first check whether we actually need to mutate this allele
    if (uniform(rng_r) < mu)
    {
        // normal distribution with as mean value the bias
        // and as sd sdmu
        std::normal_distribution <double> mutational_effect_distribution(bias, sdmu);

        G += mutational_effect_distribution(rng_r);
    }
}

// write the parameters to the DataFile
void WriteParameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< std::endl
		<< "type;" << "gonochorist_fisherian" << ";" << std::endl
		<< "popsize_init;" << N << ";" << std::endl
		<< "n_mate_sample;" << N_mate_sample << ";"<< std::endl
		<< "init_t;" << init_t << ";"<< std::endl
		<< "init_p;" << init_p << ";"<< std::endl
		<< "a;" <<  a << ";"<< std::endl
		<< "b;" <<  b << ";"<< std::endl
		<< "c;" <<  c << ";"<< std::endl
		<< "pref;" <<  pref << ";"<< std::endl
		<< "mu_p;" <<  mu_p << ";"<< std::endl
		<< "mu_t;" <<  mu_t << ";"<< std::endl
		<< "mu_std_p;" <<  sdmu_p << ";"<< std::endl
		<< "mu_std_t;" <<  sdmu_t << ";"<< std::endl
		<< "biast;" <<  biast << ";"<< std::endl
		<< "sexlimp;" <<  sexlimp << ";"<< std::endl
		<< "sexlimt;" <<  sexlimt << ";"<< std::endl
		<< "seed;" << seed << ";"<< std::endl;
} // end WriteParameters()

// initialize all the phenotypes
void Init()
{
    //++i;  prefix: it immediately adds 1
    //i++;  postfix: only adds 1 after statement has completed

	// initialize the whole populatin
	for (int i = 0; i < Nfemales; ++i)
	{
        // initialize both diploid loci
		for (int j = 0; j < 2; ++j)
		{
			Females[i].t[j] = init_t;
			Females[i].p[j] = init_p;
		}
        
        // and the expressed values
        Females[i].t_expr = init_t;
        Females[i].p_expr = init_p;
			
	}

    // initialize the male part of the population
	for (int i = 0; i < Nmales; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			Males[i].t[j] = init_t;
			Males[i].p[j] = init_p;
		}
			
        Males[i].t_expr = init_t;
        Males[i].p_expr = init_p;
	}
} // end Init

// create an offspring 
void Create_Kid(int const mother, int const father, Individual &kid)
{
	assert(mother >= 0);
    assert(mother < fsurvivors);
	assert(father >= 0); 
    assert(father < msurvivors);

    // inherit ornament
	kid.t[0] = FemaleSurvivors[mother].t[segregator(rng_r)];
    mutate(kid.t[0], mu_t, sdmu_t, -biast);

	kid.t[1] = MaleSurvivors[father].t[segregator(rng_r)];
    mutate(kid.t[1], mu_t, sdmu_t, -biast);

    // inherit preference
	kid.p[0] = FemaleSurvivors[mother].p[segregator(rng_r)];
    mutate(kid.p[0], mu_p, sdmu_p);

	kid.p[1] = MaleSurvivors[father].p[segregator(rng_r)];
    mutate(kid.p[1], mu_p, sdmu_p);
} // end Create_Kid

// survival stage
void Survive(std::ofstream &DataFile)
{
    // keep track of the 
    // number of female breeders
    fsurvivors = 0;     
    
    // store individual fitness values
    double w = 0; 

    // keep track of the 
    // mean ornament size of surviving males
    // necessary for absolute/relative preference
    // functions
    meanornsurv = 0;

    // allow females to survive
	for (int i = 0; i < Nfemales; ++i)
	{
		double p_expr = Females[i].p_expr;
		double t_expr = Females[i].t_expr;

		w = exp(-b*p_expr*p_expr + (1-sexlimt)*(-c)*t_expr*t_expr);

        // if individuals survive
        // take stats and add them to pool of survivors
        if (uniform(rng_r) < w)
        {
            FemaleSurvivors[fsurvivors++] = Females[i];
        }
	}

    msurvivors = 0;

    // male survival
	for (int i = 0; i < Nmales; ++i)
	{
		double p_expr = Males[i].p_expr;
		double t_expr = Males[i].t_expr;

		w = exp(-c*t_expr*t_expr + (1-sexlimp)*-b*p_expr*p_expr);
        
        if (uniform(rng_r) < w)
        {
            // in case of relative preferences get the mean ornament
            meanornsurv += t_expr;
            MaleSurvivors[msurvivors++] = Males[i];
        }
	}

    // extinction?
    if (fsurvivors == 0 || msurvivors == 0)
    {
        WriteParameters(DataFile);

        exit(1);
    }

    // take the average of the surviving male trait value
    meanornsurv /= msurvivors;

    assert(fsurvivors > 0);
    assert(fsurvivors < popsize);
    assert(msurvivors > 0);
    assert(msurvivors < popsize);
} // end survival stage

// mate choice

void Choose(double p, int &father) 
{
	// make arrays that hold the values of the sample of assessed males
	double Fitness[N_mate_sample];
	int Candidates[N_mate_sample];

	double sumFitness = 0;

	// check if there are enough other males to choose from
	// otherwise restrict the sample to the amount of individuals present
	int current_mate_sample = N_mate_sample > msurvivors ? msurvivors : N_mate_sample;

    std::uniform_int_distribution <int> msurvivor_sampler(0, msurvivors - 1);

    // mate choice among the sample of males
	for (int j = 0; j < current_mate_sample; ++j)
	{
		// get a random surviving male
		int random_mate = msurvivor_sampler(rng_r);

        assert(random_mate >= 0);

        assert(random_mate < msurvivors);

        // obtain a male's ornament
		double trait = MaleSurvivors[random_mate].t_expr;

        // value of the preference function
        double po = 0;

        switch(pref)
        {
            // open-ended preferences
            case 0: 
            {
                po = exp(a * trait * p);
            } break;

            // absolute preferences
            case 1:
            {
                po = exp(-a*(trait - p)*(trait - p));
            } break;

            // relative preferences
            case 2:
            {
                double reltr = trait - meanornsurv;
                po = exp(-a*(reltr - p)*(reltr - p));
            } break;

            default:
            {
                std::cout << "this should not happen, wrong preference value." << std::endl;
                exit(1);
            }
        }

        // prevent the exponential of going to infinity
        // which would make things really awkward
        if (po > 100)
        {
            po = 100;
        }

        // make a cumulative distribution of the male's
        // fitness value
		Fitness[j] = sumFitness + po;
		sumFitness=Fitness[j];

		Candidates[j] = random_mate;
	} // end for (int j = 0 

    // sample from the cumulative distribution
	double r = uniform(rng_r)*sumFitness;

    // by default mate randomly
    // e.g., if the cumulative distribution is flat
	father = msurvivor_sampler(rng_r);

    // probability that a male is chosen is proportional
    // to the size of his share in the cumulative distribution
	for (int j = 0; j < current_mate_sample; ++j)
	{
		if (r <= Fitness[j])
		{
			father=Candidates[j];				
			break;	
		}
	}

    assert(father >= 0 && father < msurvivors);

} // end ChooseMates
// produce the next generation
void NextGen()
{
    int offspring = 0;

    // do some stats of the adult reproductive
    // population
    if (do_stats)
    {
        for (int i = 0; i < msurvivors; ++i)
        {
            father_eggs[i] = 0;
        }

        for (int i = 0; i < fsurvivors; ++i)
        {
            mother_eggs[i] = 0;
        }
    } // end if(do_stats)

    // let the surviving females choose a mate
	for (int i = 0; i < fsurvivors; ++i)
	{
		int Father = -1;
        
        double expr_p = FemaleSurvivors[i].p_expr;

		Choose(expr_p, Father);

		assert(Father >= 0 && Father < msurvivors);

        // do some stats on the mating pairs
        if (do_stats)
        {
            mother_eggs[i] += clutch_size;
            father_eggs[Father] += clutch_size;
        }

        // for each offspring to be produced
        // store the indices of the parents
        // we then make offspring later
        for (int j = 0; j < clutch_size; ++j)
        {
            Parents[offspring][0] = i;
            Parents[offspring][1] = Father;
            ++offspring;
        }
	} // end for (int i = 0; i < fsurvivors

    int sons = 0;
    int daughters = 0;

    // the size of the resulting population in the next 
    // generation is dependent on the number of 
    // offspring produced
    popsize = offspring < N ? offspring : N;

   std::uniform_int_distribution <int> offspring_sampler(0, offspring - 1);

    // replace the next generation
    for (int i = 0; i < popsize; ++i)
    {
        // create an offspring
        Individual Kid;

        int random_offspring = offspring_sampler(rng_r);

        // randomly sample an offspring to replace the population
        Create_Kid(Parents[random_offspring][0], Parents[random_offspring][1], Kid);
            
        double t = 0.5 * ( Kid.t[0] + Kid.t[1]);
        double p = 0.5 * ( Kid.p[0] + Kid.p[1]);

        Kid.t_expr = t; 
        Kid.p_expr = p; 

        assert(random_offspring >= 0);
        assert(random_offspring < N*clutch_size);


        assert(Parents[random_offspring][0] >= 0); 
        assert(Parents[random_offspring][0] < fsurvivors);
        assert(Parents[random_offspring][1] >= 0);
        assert(Parents[random_offspring][1] < msurvivors);

        // it's a boy
        if (uniform(rng_r) < 0.5)
        {
            Males[sons++] = Kid;
        }
        else
        {
            Females[daughters++] = Kid;
        }
    }

    Nmales = sons;
    Nfemales = daughters;
}



// write the data
void WriteData(std::ofstream &DataFile)
{
    // in case population is extinct,
    // quit
	if (Nmales == 0 || Nfemales == 0)
	{
		WriteParameters(DataFile);
		exit(1);
	}

    double meanp = 0.0;
    double meant = 0.0; 
    double ssp = 0.0;
    double sst = 0.0;
    double spt = 0.0;

	double p,t,meanmrs,meanfrs,varfrs,varmrs;
	double ssmrs = 0, ssfrs = 0, summrs=0, sumfrs=0;

    // calculate means and variances for the males
	for (int i = 0; i < Nmales; ++i)
	{
		p = 0.5 * ( Males[i].p[0] + Males[i].p[1]);
		t = 0.5 * ( Males[i].t[0] + Males[i].t[1]);

        meanp += p;
        meant += t;

        ssp += p * p;
        sst += t * t;
        spt += t * p;

        if (i < msurvivors)
        {
            summrs += father_eggs[i];
            ssmrs += father_eggs[i] * father_eggs[i];
        }
	}

    // calculate means and variances for the females
	for (int i = 0; i < Nfemales; ++i)
	{
		p = 0.5 * ( Females[i].p[0] + Females[i].p[1]);
		t = 0.5 * ( Females[i].t[0] + Females[i].t[1]);

        meanp += p;
        meant += t;

        ssp += p * p;
        sst += t * t;
        spt += t * p;
        
        if (i < fsurvivors)
        {
            sumfrs += mother_eggs[i];
            ssfrs += mother_eggs[i] * mother_eggs[i];
        }
	} 

    double varp = 0; 
    double vart = 0; 
    double covpt = 0;

    double sum_sexes = Nmales + Nfemales;

    meant /= sum_sexes;
    meanp /= sum_sexes;
    vart = sst / sum_sexes - meant * meant;
    varp = ssp / sum_sexes - meanp * meanp;
    covpt = spt / sum_sexes - meanp * meant;

    meanfrs = sumfrs / Nfemales;
    varfrs = ssfrs / Nfemales - meanfrs * meanfrs;

    meanmrs = summrs / Nmales;
    varmrs = ssmrs / Nmales - meanmrs * meanmrs;


    // output of all the statistics
	DataFile << generation 
		<< ";" << meanp
		<< ";" << meant
		
        << ";" << meanfrs 
        << ";" << meanmrs 
		
        << ";" << varp
		<< ";" << vart
		<< ";" << covpt 
		
        << ";" << varfrs 
        << ";" << varmrs 

        << ";" << msurvivors
        << ";" << fsurvivors
		<< ";" << sum_sexes
	       	<< ";" << std::endl;
} // end WriteData

// headers of the datafile
void WriteDataHeaders(std::ofstream &DataFile)
{
	DataFile << "generation" 
		<< ";mean_p" 
		<< ";mean_t" 
        << ";meanfrs"
        << ";meanmrs"
        << ";varp"
        << ";vart"
        << ";covpt"
        << ";varfrs"
        << ";varmrs"
        << ";surviving_males"
        << ";surviving_females"
		<< ";N;"
		<< std::endl;
} // end WriteDataHeaders

// the core part of the code
int main(int argc, char ** argv)
{
    // initialize parameters based on command line arguments
	initArguments(argc, argv);

    // initialize output file
    std::ofstream output_file(file_name.c_str());

    // write the headers to the data file
	WriteDataHeaders(output_file);

    // initialize the population
	Init();

    // loop through each generation
	for (generation = 0; generation <= NumGen; ++generation)
	{
        // decide whether we want to write out stats or not
		do_stats = generation % skip == 0;

        // individuals survive (or not)
		Survive(output_file);
        
        // individuals reproduce (or not)
        // female choice, so all females reproduce (by definition)
        // but not all males.
        NextGen();
        
        if (do_stats)
		{
			WriteData(output_file);
		}
	}

	WriteParameters(output_file);
} // end main
