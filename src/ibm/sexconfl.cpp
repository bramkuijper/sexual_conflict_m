// sexual conflict with maternal and paternal effects
// 
//
//     Copyright (C) 2017 Bram Kuijper
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


// load some libraries 
#define NDEBUG
//#define DISTRIBUTION
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>

// libraries for random number generations
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

using namespace std;

// define all the parameters
const int N = 5000; // population size
const int N_mate_sample = 10; // how many mates does each female encounter?
const int clutch_size = 20; // how many eggs does she produce?
double a = 1.0; // strength of selection against non-optimal # matings
double cs = 0.5; // cost of sensitivity
double ct = 0.5; // cost of threshold
double co = 0.2; // cost of offense trait
double const thr_opt = 2; // optimal values for each of the traits
double const sen_opt = 3;
double const off_opt = 2;
double const init_sen = sen_opt; // specify initial values for each trait
double const init_off = off_opt; 
double const init_thr = thr_opt; 
double gammathr = 2; // strengths of selection
double gammasen = 2;
double const powoff = 2;
double mu_off 	  = 0.05;            // mutation rate
double mu_off_p     = 0.05;            // mutation rate
double mu_thr 	  = 0.05;            // mutation rate
double mu_thr_m     = 0.05;            // mutation rate
double mu_sen     = 0.05;            // mutation rate
double mu_sen_m     = 0.05;            // mutation rate
double sdmu_off         = 0.4;			 // standard deviation mutation size
double sdmu_off_p         = 0.4;			 // standard deviation mutation size
double sdmu_thr         = 0.4;			 // standard deviation mutation size
double sdmu_thr_m         = 0.4;			 // standard deviation mutation size
double sdmu_sen         = 0.4;			 // standard deviation mutation size
double sdmu_sen_m         = 0.4;			 // standard deviation mutation size
const double NumGen = 50000;
const int skip = 10;
double theta_psi = 1;

bool do_stats = 0;

int generation = 0;
int Nfemales = N / 2, Nmales = N / 2;
int popsize = 0;
unsigned seed = 0;
int msurvivors = 0;
int fsurvivors = 0;
int unmated_f = 0;
int nmatings_per_female = 0;

int father_eggs[N];
int mother_eggs[N];

int mpartners[N];
int fpartners[N];

// the individual struct
struct Individual
{
	double off[2]; // male offense trait
	double thr[2]; // female threshold
	double sen[2]; // female sensitivity
	
    double off_p[2]; // male offense trait paternal effect
	double thr_m[2]; // female threshold maternal effect
	double sen_m[2]; // female sensitivity maternal effect

    double e_off;
    double e_thr;
    double e_sen;
};

typedef Individual Population[N];
Population Females, Males, FemaleSurvivors, MaleSurvivors;
typedef Individual Kids[N*clutch_size];
Kids NewPop;

string filename("sim_sexualconflict");
string filename_new(create_filename(filename));
ofstream DataFile(filename_new.c_str());  // output file 

#ifdef DISTRIBUTION
string filename_new2(create_filename("sim_sexualconflict"));
ofstream distfile(filename_new2.c_str());
#endif //DISTRIBUTION

void initArguments(int argc, char *argv[])
{
	a = atof(argv[1]);
	co = atof(argv[2]);
	ct = atof(argv[3]);
	cs = atof(argv[4]);
	mu_off = atof(argv[5]);
	mu_thr = atof(argv[6]);
	mu_sen = atof(argv[7]);
	mu_off_p = atof(argv[8]);
	mu_thr_m = atof(argv[9]);
	mu_sen_m = atof(argv[10]);
	sdmu_off = atof(argv[11]);
	sdmu_thr = atof(argv[12]);
	sdmu_sen = atof(argv[13]);
	sdmu_off_p = atof(argv[14]);
	sdmu_thr_m = atof(argv[15]);
	sdmu_sen_m = atof(argv[16]);
    theta_psi = atof(argv[17]);
}

void Mutate(double &G, double const mu, double const sdmu)
{
	G += gsl_rng_uniform(r)<mu ? gsl_ran_gaussian(r, sdmu) : 0;
}


void WriteParameters()
{
	DataFile << endl
		<< endl
		<< "type;" << "sexual_conflict" << ";" << endl
		<< "popsize_init;" << N << ";" << endl
		<< "n_mate_sample;" << N_mate_sample << ";"<< endl
		<< "a;" <<  a << ";"<< endl
		<< "coff;" <<  co << ";"<< endl
		<< "cthr;" <<  ct << ";"<< endl
		<< "csen;" <<  cs << ";"<< endl
		<< "thr_opt;" <<  thr_opt << ";"<< endl
		<< "sen_opt;" <<  sen_opt << ";"<< endl
		<< "off_opt;" <<  off_opt << ";"<< endl
		<< "theta_psi;" << theta_psi << ";"<< endl
		<< "mu_off;" <<  mu_off << ";"<< endl
		<< "mu_thr;" <<  mu_thr << ";"<< endl
		<< "mu_sen;" <<  mu_sen << ";"<< endl
		<< "mu_off_p;" <<  mu_off_p << ";"<< endl
		<< "mu_thr_m;" <<  mu_thr_m << ";"<< endl
		<< "mu_sen_m;" <<  mu_sen_m << ";"<< endl
		<< "mu_std_off;" <<  sdmu_off << ";"<< endl
		<< "mu_std_thr;" <<  sdmu_thr << ";"<< endl
		<< "mu_std_sen;" <<  sdmu_sen << ";"<< endl
		<< "mu_std_off_p;" <<  sdmu_off_p << ";"<< endl
		<< "mu_std_thr_m;" <<  sdmu_thr_m << ";"<< endl
		<< "mu_std_sen_m;" <<  sdmu_sen_m << ";"<< endl
		<< "seed;" << seed << ";"<< endl;
}

// initialize all the phenotypes
void Init()
{
	seed = get_nanoseconds();

    // set up the random number generators
    // (from the gnu gsl library)
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    
    for (int i = 0; i < N/2; ++i)
    {
        Females[i].off[0] = init_off;
        Females[i].off[1] = init_off;
        Females[i].sen[0] = init_sen;
        Females[i].sen[1] = init_sen;
        Females[i].thr[0] = init_thr;
        Females[i].thr[1] = init_thr;
        
        Males[i].off[0] = init_off;
        Males[i].off[1] = init_off;
        Males[i].sen[0] = init_sen;
        Males[i].sen[1] = init_sen;
        Males[i].thr[0] = init_thr;
        Males[i].thr[1] = init_thr;
    
        Females[i].e_thr =  init_thr;
        Females[i].e_sen =  init_sen;
        Males[i].e_off =  init_off;
    }

    Nmales = N/2;
    Nfemales = N/2; 

    popsize = N;
}

void Create_Kid(int mother, int father, Individual &kid)
{
	assert(mother >= 0 && mother < fsurvivors);
	assert(father >= 0 && father < msurvivors);

    // inherit offense trait
    kid.off[0] = FemaleSurvivors[mother].off[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.off[0], mu_off, sdmu_off);
    kid.off[1] = MaleSurvivors[father].off[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.off[1], mu_off, sdmu_off);
    
    // inherit offense trait paternal effect genes
    kid.off_p[0] = FemaleSurvivors[mother].off_p[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.off_p[0], mu_off_p, sdmu_off_p);
    kid.off_p[1] = MaleSurvivors[father].off_p[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.off_p[1], mu_off_p, sdmu_off_p);

    // inherit threshold
    kid.thr[0] = FemaleSurvivors[mother].thr[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.thr[0], mu_thr, sdmu_thr);
    kid.thr[1] = MaleSurvivors[father].thr[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.thr[1], mu_thr, sdmu_thr);
    
    // inherit threshold maternal effect genes
    kid.thr_m[0] = FemaleSurvivors[mother].thr_m[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.thr_m[0], mu_thr_m, sdmu_thr_m);
    kid.thr_m[1] = MaleSurvivors[father].thr_m[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.thr_m[1], mu_thr_m, sdmu_thr_m);

    // inherit sensitivity
    kid.sen[0] = FemaleSurvivors[mother].sen[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.sen[0], mu_sen, sdmu_sen);
    kid.sen[1] = MaleSurvivors[father].sen[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.sen[1], mu_sen, sdmu_sen);

    // inherit sen_msitivity maternal effect genes
    kid.sen_m[0] = FemaleSurvivors[mother].sen_m[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.sen_m[0], mu_sen_m, sdmu_sen_m);
    kid.sen_m[1] = MaleSurvivors[father].sen_m[gsl_rng_uniform_int(r, 2)];
    Mutate(kid.sen_m[1], mu_sen_m, sdmu_sen_m);

    // express phenotypes
    kid.e_off = 0.5 * (kid.off[0] + kid.off[1]) 
        + 0.5 * (kid.off_p[0] + kid.off_p[1]) * MaleSurvivors[father].e_off;
    
    kid.e_thr = 0.5 * (kid.thr[0] + kid.thr[1]) 
        + 0.5 * (kid.thr_m[0] + kid.thr_m[1]) * FemaleSurvivors[mother].e_thr;
    
    kid.e_sen = 0.5 * (kid.sen[0] + kid.sen[1]) 
        + 0.5 * (kid.sen_m[0] + kid.sen_m[1]) * FemaleSurvivors[mother].e_sen;
}


// threshold and sensitivity affect survival
// psi affects fecundity
// we vary this later
void Survive()
{
    fsurvivors = 0;
    double w = 0;

    // let females survive
	for (int i = 0; i < Nfemales; ++i)
	{
        // we use the same scaling of fitness as Rowe et al 
        // but we use exponents here:
        w = exp(-ct * pow(Females[i].e_thr - thr_opt,2) 
                    - cs * pow(Females[i].e_sen - sen_opt,2));

        // female is surviving
        if (gsl_rng_uniform(r) < w)
        {
            FemaleSurvivors[fsurvivors++] = Females[i];
        }
	}

    msurvivors = 0;

	// now get the male trait values into a viability measure
	for (int i = 0; i < Nmales; ++i)
	{
        // natural selection on a male offense trait
		w = exp(-co * pow(Males[i].e_off - off_opt, 2));
        
        // survival of males
        if (gsl_rng_uniform(r) < w)
        {
            MaleSurvivors[msurvivors++] = Males[i];
        }
	}

    if (fsurvivors == 0 || msurvivors == 0)
    {
        WriteParameters();

        exit(1);
    }

    assert(fsurvivors > 0 && fsurvivors < popsize);
    assert(msurvivors > 0 && msurvivors < popsize);

}

void Choose(int &mother, int &offspring) 
{
	// check if there are enough other individuals to choose from
	// otherwise restrict the sample to the amount of individuals present
	int encounters = N_mate_sample > msurvivors ? msurvivors : N_mate_sample;

    int nMatings = 0;

    int mates[N_mate_sample];

    // encounter a number of males 
    // mate with them only if the threshold and sensitivity allow it
    // when they mate this will be costly
	for (int i = 0; i < encounters; ++i)
	{
		// get a random male survivor
		int random_mate = gsl_rng_uniform_int(r, msurvivors);
        
        // calculate the mating rate
        // see p. 57 2nd column in Rowe et al
        double signal = FemaleSurvivors[mother].e_sen * (
                        MaleSurvivors[random_mate].e_off 
                        - FemaleSurvivors[mother].e_thr
                        );

//        double psi = 1.0 / (1.0 + exp(-signal)); // no signal: psi = 0.5
//
//      In Rowe's model theta_psi is always .5. Here I account for the fact
//      that other optima may be possible as well.
        double psi = 1.0 / (1.0 + ((1.0/theta_psi) - 1.0) * exp(-signal)); // no signal: psi = 0.5

        // will this male mate with this female?
        if (gsl_rng_uniform(r) < psi)
        {
            mates[nMatings++]=random_mate;
        }
	}

    // the more the total number of matings diverges from the optimum, the lower
    // the fecundity
    double fecundity_continuous = (double) clutch_size * 
        exp(-a * pow((double) nMatings / encounters - theta_psi, 2.0));

    int fecundity = floor(fecundity_continuous);

    // rounding of a continuous fecundity number to a discrete one
    double delta_round = fecundity_continuous - fecundity;

    if (gsl_rng_uniform(r) < delta_round)
    {
        ++fecundity;
    }

    assert(fecundity <= clutch_size);

    // ok, now make offspring where paternity is 
    // randomly distributed over mates
    if (nMatings > 0)
    {
        for (int i = 0; i < fecundity; ++i)
        {
            int randnum = gsl_rng_uniform_int(r, nMatings);
            int current_father = mates[randnum];

            Individual Kid;
            Create_Kid(mother, current_father, Kid);
            NewPop[offspring++] = Kid;
        }
    }
    else
    {
        ++unmated_f; // keep track of unmated females
    }
} // end ChooseMates

// choose the mate
// this is a simple unilateral mate choice function
// later on we have to think about mutual mate choice
void NextGen()
{
    int offspring = 0;

    unmated_f = 0;
    
    // let the surviving females choose a mate
	for (int i = 0; i < fsurvivors; ++i)
	{
		Choose(i, offspring);
	}
    
    if (offspring == 0)
    {
        WriteParameters();
        exit(1);
    }

    assert(offspring > 0 && offspring < popsize*clutch_size);

    int sons = 0;
    int daughters = 0;

    popsize = offspring < N ? offspring : N;

    // replace population with offspring
    for (int i = 0; i < popsize; ++i)
    {
        if (gsl_rng_uniform(r) < 0.5)
        {
            Males[sons++] = NewPop[gsl_rng_uniform_int(r, offspring)];
        }
        else
        {
            Females[daughters++] = NewPop[gsl_rng_uniform_int(r, offspring)];
        }
    }

    Nmales = sons;
    Nfemales = daughters;
}

void WriteData()
{
	if (Nmales == 0 || Nfemales == 0)
	{
		WriteParameters();
		exit(1);
	}

    double off, sen, thr, off_p, sen_m, thr_m;

    double mean_e_off = 0;
    double meanoff = 0;
    double mean_e_sen = 0;
    double meansen = 0;
    double mean_e_thr = 0;
    double meanthr = 0;

    double ss_e_off = 0;
    double ssoff = 0;
    double ss_e_sen = 0;
    double sssen = 0;
    double ss_e_thr = 0;
    double ssthr = 0;

    double meanoff_p = 0;
    double meansen_m = 0;
    double meanthr_m = 0;

    double ssoff_p = 0;
    double sssen_m = 0;
    double ssthr_m = 0;

	for (int i = 0; i < Nmales; ++i)
	{
        off = 0.5 * ( Males[i].off[0] + Males[i].off[1]);
        sen = 0.5 * ( Males[i].sen[0] + Males[i].sen[1]);
        thr = 0.5 * ( Males[i].thr[0] + Males[i].thr[1]);
        
        off_p = 0.5 * ( Males[i].off_p[0] + Males[i].off_p[1]);
        sen_m = 0.5 * ( Males[i].sen_m[0] + Males[i].sen_m[1]);
        thr_m = 0.5 * ( Males[i].thr_m[0] + Males[i].thr_m[1]);

        meanoff += off;
        meansen += sen;
        meanthr += thr;

        mean_e_off += Males[i].e_off;
        mean_e_sen += Males[i].e_sen;
        mean_e_thr += Males[i].e_thr;
        
        meanoff_p += off_p;
        meansen_m += sen_m;
        meanthr_m += thr_m;
        
        ssoff += off*off;
        sssen += sen*sen;
        ssthr += thr*thr;
        
        ssoff_p += off_p*off_p;
        sssen_m += sen_m*sen_m;
        ssthr_m += thr_m*thr_m;

        ss_e_off += Males[i].e_off * Males[i].e_off;
        ss_e_sen += Males[i].e_sen * Males[i].e_sen;
        ss_e_thr += Males[i].e_thr * Males[i].e_thr;
	}

	for (int i = 0; i < Nfemales; ++i)
	{
        off = 0.5 * (Females[i].off[0] + Females[i].off[1]);
        sen = 0.5 * (Females[i].sen[0] + Females[i].sen[1]);
        thr = 0.5 * (Females[i].thr[0] + Females[i].thr[1]);
        
        off_p = 0.5 * ( Females[i].off_p[0] + Females[i].off_p[1]);
        sen_m = 0.5 * ( Females[i].sen_m[0] + Females[i].sen_m[1]);
        thr_m = 0.5 * ( Females[i].thr_m[0] + Females[i].thr_m[1]);

        meanoff += off;
        meansen += sen;
        meanthr += thr;
        
        mean_e_off += Females[i].e_off;
        mean_e_sen += Females[i].e_sen;
        mean_e_thr += Females[i].e_thr;
        
        meanoff_p += off_p;
        meansen_m += sen_m;
        meanthr_m += thr_m;
        
        ssoff += off*off;
        sssen += sen*sen;
        ssthr += thr*thr;
        
        ssoff_p += off_p*off_p;
        sssen_m += sen_m*sen_m;
        ssthr_m += thr_m*thr_m;
        
        ss_e_off += Females[i].e_off * Females[i].e_off;
        ss_e_sen += Females[i].e_sen * Females[i].e_sen;
        ss_e_thr += Females[i].e_thr * Females[i].e_thr;
	}

    double sum_sexes = Nmales + Nfemales;

    meanoff /= sum_sexes;
    meansen /= sum_sexes;
    meanthr /= sum_sexes;
    mean_e_off /= sum_sexes;
    mean_e_sen /= sum_sexes;
    mean_e_thr /= sum_sexes;
    meanoff_p /= sum_sexes;
    meansen_m /= sum_sexes;
    meanthr_m /= sum_sexes;

	DataFile << generation << ";"
            << meanoff << ";" 
            << meansen << ";" 
            << meanthr << ";" 
            << mean_e_off << ";" 
            << mean_e_sen << ";" 
            << mean_e_thr << ";" 
            << meanoff_p << ";" 
            << meansen_m << ";" 
            << meanthr_m << ";" 
            << ssoff / sum_sexes - meanoff * meanoff << ";" 
            << sssen / sum_sexes  - meansen * meansen << ";" 
            << ssthr / sum_sexes - meanthr * meanthr << ";" 
            << ss_e_off / sum_sexes - mean_e_off * mean_e_off << ";" 
            << ss_e_sen / sum_sexes  - mean_e_sen * mean_e_sen << ";" 
            << ss_e_thr / sum_sexes - mean_e_thr * mean_e_thr << ";" 
            << ssoff_p / sum_sexes - meanoff_p * meanoff_p << ";" 
            << sssen_m / sum_sexes  - meansen_m * meansen_m << ";" 
            << ssthr_m / sum_sexes - meanthr_m * meanthr_m << ";" 
            << Nmales << ";"
            << Nfemales << ";"
            << msurvivors << ";"
            << fsurvivors << ";"
            << unmated_f << ";"
            << endl;

}

void WriteDataHeaders()
{
	DataFile << "generation;" << 
        "meanoff;" << 
        "meansen;" << 
        "meanthr;" << 
        "meanoff_phen;" <<
        "meansen_phen;" << 
        "meanthr_phen;" << 
        "meanoff_p;" <<
        "meansen_m;" << 
        "meanthr_m;" << 
        "varoff;" << 
        "varsen;" << 
        "varthr;" << 
        "varoff_phen;" << 
        "varsen_phen;" << 
        "varthr_phen;" << 
        "varoff_p;" << 
        "varsen_m;" << 
        "varthr_m;" << 
        "Nm;" <<
        "Nf;" <<
        "Nmsurv;" <<
        "Nfsurv;" <<
        "Nfunmated;" << endl;

}

int main(int argc, char ** argv)
{
	initArguments(argc, argv);
	WriteDataHeaders();
	Init();

	for (generation = 0; generation <= NumGen; ++generation)
	{
		do_stats = generation % skip == 0;

		Survive();
		
        NextGen();
        
        if (do_stats)
		{
			WriteData();
		}
	}

	WriteParameters();
}
