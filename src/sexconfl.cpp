// sexual conflict with multiple offense and resistance traits
// 
//
//     Copyright (C) 2010 Bram Kuijper
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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "bramauxiliary.h"


// random number generator 
// see http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Generation.html#Random-Number-Generation 
gsl_rng_type const * T; // gnu scientific library rng type
gsl_rng *r; // gnu scientific rng 

using namespace std;

const int N = 5000;
const int N_mate_sample = 10;
const int clutch_size = 20;
double a = 1.0; // strength of selection against non-optimal # matings
double cs = 0.5; // cost of sensitivity
double ct = 0.5; // cost of     
double co = 0.2; // cost of trait
double const thr_opt = 2;
double const sen_opt = 3;
double const off_opt = 2;
double const init_sen = sen_opt; // traits start at their optima
double const init_off = off_opt; 
double const init_thr = thr_opt; 
double ipowthr = 0.1;
double ipowsen = 0.1;
double gammathr = 2;
double gammasen = 2;
double const powoff = 2;
double mu_off 	  = 0.05;            // mutation rate
double mu_thr 	  = 0.05;            // mutation rate
double mu_sen     = 0.05;            // mutation rate
double sdmu_off         = 0.4;			 // standard deviation mutation size
double sdmu_thr         = 0.4;			 // standard deviation mutation size
double sdmu_sen         = 0.4;			 // standard deviation mutation size
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

    // TODO
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
	sdmu_off = atof(argv[8]);
	sdmu_thr = atof(argv[9]);
	sdmu_sen = atof(argv[10]);
    theta_psi = atof(argv[11]);
    ipowthr = atof(argv[12]);
    ipowsen = atof(argv[13]);
}

void MutateOff(double &G)
{
	G += gsl_rng_uniform(r)<mu_off ? gsl_ran_gaussian(r, sdmu_off) : 0;
}

void MutateThr(double &G)
{
	G+= gsl_rng_uniform(r)<mu_thr ? gsl_ran_gaussian(r, sdmu_thr) : 0; 
}

void MutateSen(double &G)
{
	G+= gsl_rng_uniform(r)<mu_sen ? gsl_ran_gaussian(r, sdmu_sen) : 0; 
}

void WriteParameters()
{
	DataFile << endl
		<< endl
		<< "type:;" << "sexual_conflict" << ";" << endl
		<< "popsize_init:;" << N << ";" << endl
		<< "n_mate_sample:;" << N_mate_sample << ";"<< endl
		<< "a:;" <<  a << ";"<< endl
		<< "coff:;" <<  co << ";"<< endl
		<< "cthr:;" <<  ct << ";"<< endl
		<< "csen:;" <<  cs << ";"<< endl
		<< "thr_opt:;" <<  thr_opt << ";"<< endl
		<< "sen_opt:;" <<  sen_opt << ";"<< endl
		<< "off_opt:;" <<  off_opt << ";"<< endl
		<< "powthr:;" <<  ipowthr << ";"<< endl
		<< "powsen:;" <<  ipowsen << ";"<< endl
		<< "theta_psi:;" << theta_psi << ";"<< endl
		<< "mu_off:;" <<  mu_off << ";"<< endl
		<< "mu_thr:;" <<  mu_thr << ";"<< endl
		<< "mu_sen:;" <<  mu_sen << ";"<< endl
		<< "mu_std_off:;" <<  sdmu_off << ";"<< endl
		<< "mu_std_thr:;" <<  sdmu_thr << ";"<< endl
		<< "mu_std_sen:;" <<  sdmu_sen << ";"<< endl
		<< "seed:;" << seed << ";"<< endl;
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
    kid.off = FemaleSurvivors[mother].off[gsl_rng_uniform_int(r, 2)];
    MutateOff(kid.off);
    kid.off = MaleSurvivors[father].off[gsl_rng_uniform_int(r, 2)];
    MutateOff(kid.off);

    // inherit threshold
    kid.thr = FemaleSurvivors[mother].thr[gsl_rng_uniform_int(r, 2)];
    MutateThr(kid.thr);
    kid.thr = MaleSurvivors[father].thr[gsl_rng_uniform_int(r, 2)];
    MutateThr(kid.thr);

    // inherit sensitivity
    kid.sen = FemaleSurvivors[mother].sen[gsl_rng_uniform_int(r, 2)];
    MutateSen(kid.sen);
    kid.sen = MaleSurvivors[father].sen[gsl_rng_uniform_int(r, 2)];
    MutateSen(kid.sen);
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
    double fecundity_continuous = (double) clutch_size * exp(-a * pow((double) nMatings / encounters - theta_psi, 2.0));

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
    
    if (do_stats)
    {
        for (int i = 0; i < msurvivors; ++i)
        {
            father_eggs[i] = 0;
            mpartners[i] = 0;
        }
        
        for (int i = 0; i < fsurvivors; ++i)
        {
            mother_eggs[i] = 0;
            fpartners[i] = 0;
        }
        
        unmated_f = 0;
    }

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

    for (int i = 0; i < popsize; ++i)
    {
        if (Uniform() < 0.5)
        {
            Males[sons] = NewPop[RandomNumber(offspring)];
   
            for (int j = 0; j < nTraits; ++j)
            {
                double off = 0.5 * ( Males[sons].off[j][0] + Males[sons].off[j][1]);
                Males[sons].e_off[j] = off; 
            }
            ++sons;
        }
        else
        {
            Females[daughters] = NewPop[RandomNumber(offspring)];
            
            for (int j = 0; j < nTraits; ++j)
            {
                double thr = 0.5 * ( Females[daughters].thr[j][0] + Females[daughters].thr[j][1]);
                Females[daughters].e_thr[j] = thr;
                double sen = 0.5 * ( Females[daughters].sen[j][0] + Females[daughters].sen[j][1]);
                Females[daughters].e_sen[j] = sen;
            }
            ++daughters;
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

    for (int i = 0; i < nTraits; ++i)
    {
        stat_reset(stat_start_phen_off[i]);
        stat_reset(stat_start_phen_thr[i]);
        stat_reset(stat_start_phen_sen[i]);
        stat_reset(stat_start_off[i]);
        stat_reset(stat_start_thr[i]);
        stat_reset(stat_start_sen[i]);
        jstat_reset(jstat_start_offsen[i]);
        jstat_reset(jstat_start_offthr[i]);
        jstat_reset(jstat_start_thrsen[i]);
    }

    int sumfrs = 0; 
    int summrs = 0; 
    int ssfrs = 0; 
    int ssmrs = 0; 

    double meanmrs, meanfrs, varfrs, varmrs;

    double f_matings = 0;
    double ss_f_matings = 0;
    double m_matings = 0;
    double ss_m_matings = 0;

	for (int i = 0; i < Nmales; ++i)
	{
        double off, sen, thr;

        for (int j = 0; j < nTraits; ++j)
        {
            off = 0.5 * ( Males[i].off[j][0] + Males[i].off[j][1]);
            sen = 0.5 * ( Males[i].sen[j][0] + Males[i].sen[j][1]);
            thr = 0.5 * ( Males[i].thr[j][0] + Males[i].thr[j][1]);

#ifdef DISTRIBUTION 
            distfile << generation << ";"
                    << setprecision(5) << off << ";"
                    << setprecision(5) << sen << ";"
                    << setprecision(5) << 0 << ";" 
                    << setprecision(5) << thr << ";" << endl;
#endif // DISTRIBUTION

            stat_addval(stat_start_off[j], off);
            stat_addval(stat_start_thr[j], thr);
            stat_addval(stat_start_sen[j], sen);
            stat_addval(stat_start_phen_off[j], Males[i].e_off[j]);
            jstat_addval(jstat_start_offsen[j], off, sen); 
            jstat_addval(jstat_start_offthr[j], off, thr); 
            jstat_addval(jstat_start_thrsen[j], sen, thr); 
        }
            
        if (i < msurvivors)
        {
            summrs += father_eggs[i];
            ssmrs += father_eggs[i] * father_eggs[i];

            m_matings += mpartners[i];
            ss_m_matings += mpartners[i] * mpartners[i];

        }
	}

	for (int i = 0; i < Nfemales; ++i)
	{
        double off, sen, thr;

        for (int j = 0; j < nTraits; ++j)
        {
            off = 0.5 * ( Females[i].off[j][0] + Females[i].off[j][1]);
            sen = 0.5 * ( Females[i].sen[j][0] + Females[i].sen[j][1]);
            thr = 0.5 * ( Females[i].thr[j][0] + Females[i].thr[j][1]);


#ifdef DISTRIBUTION 
            distfile << generation << ";"
                    << setprecision(5) << off << ";"
                    << setprecision(5) << sen << ";"
                    << setprecision(5) << 0 << ";" 
                    << setprecision(5) << thr << ";" << endl;
#endif // DISTRIBUTION
            
            stat_addval(stat_start_phen_thr[j], Females[i].e_thr[j]);
            stat_addval(stat_start_phen_sen[j], Females[i].e_sen[j]);
            stat_addval(stat_start_off[j], off);
            stat_addval(stat_start_thr[j], thr);
            stat_addval(stat_start_sen[j], sen);
            jstat_addval(jstat_start_offsen[j], off, sen); 
            jstat_addval(jstat_start_offthr[j], off, thr); 
            jstat_addval(jstat_start_thrsen[j], sen, thr); 
        }

        if (i < fsurvivors)
        {
            sumfrs += mother_eggs[i];
            ssfrs += mother_eggs[i] * mother_eggs[i];
            f_matings += fpartners[i];
            ss_f_matings += fpartners[i] * fpartners[i];
        }
	}

    for (int i = 0; i < nTraits; ++i)
    {
        stat_finalize(stat_start_phen_thr[i]);
        stat_finalize(stat_start_thr[i]);
        stat_finalize(stat_ns_thr[i]);

        stat_finalize(stat_start_phen_sen[i]);
        stat_finalize(stat_start_sen[i]);
        stat_finalize(stat_ns_sen[i]);
        
        stat_finalize(stat_start_phen_off[i]);
        stat_finalize(stat_start_off[i]);
        stat_finalize(stat_ns_off[i]);

        jstat_finalize(jstat_start_offsen[i], stat_start_off[i].mean, stat_start_sen[i].mean);
        jstat_finalize(jstat_start_offthr[i], stat_start_off[i].mean, stat_start_thr[i].mean);
        jstat_finalize(jstat_start_thrsen[i], stat_start_sen[i].mean, stat_start_thr[i].mean);
    }
    
    double sum_sexes = Nmales + Nfemales;

    meanfrs = (double) sumfrs / Nfemales;
    meanmrs = (double) summrs / Nmales;
    varfrs = (double) ssfrs / Nfemales - meanfrs * meanfrs;
    varmrs = (double) ssmrs / Nmales - meanmrs * meanmrs;

    m_matings /= (double) msurvivors;
    f_matings /= (double) fsurvivors;

    double var_m_matings = ss_m_matings / msurvivors - m_matings * m_matings;
    double var_f_matings = ss_f_matings / fsurvivors - f_matings * f_matings;

	DataFile << generation;

    for (int i = 0; i < nTraits; ++i)
    {
		DataFile << ";" << setprecision(5) << stat_start_off[i].mean
		<< ";" << setprecision(5) << stat_start_thr[i].mean
		<< ";" << setprecision(5) << stat_start_sen[i].mean
		
		<< ";" << setprecision(5) << stat_start_off[i].mean_ci
		<< ";" << setprecision(5) << stat_start_thr[i].mean_ci
		<< ";" << setprecision(5) << stat_start_sen[i].mean_ci

		<< ";" << setprecision(5) << stat_start_phen_off[i].mean
		<< ";" << setprecision(5) << stat_start_phen_thr[i].mean
		<< ";" << setprecision(5) << stat_start_phen_sen[i].mean
		
		<< ";" << setprecision(5) << stat_ns_off[i].mean
		<< ";" << setprecision(5) << stat_ns_thr[i].mean
		<< ";" << setprecision(5) << stat_ns_sen[i].mean
		
        << ";" << setprecision(5) << stat_ns_off[i].mean_ci
		<< ";" << setprecision(5) << stat_ns_thr[i].mean_ci
		<< ";" << setprecision(5) << stat_ns_sen[i].mean_ci

		<< ";" << setprecision(5) << stat_start_off[i].var
		<< ";" << setprecision(5) << stat_start_thr[i].var
		<< ";" << setprecision(5) << stat_start_sen[i].var
		
        << ";" << setprecision(5) << stat_start_phen_off[i].var
		<< ";" << setprecision(5) << stat_start_phen_thr[i].var
		<< ";" << setprecision(5) << stat_start_phen_sen[i].var

		<< ";" << setprecision(5) << stat_ns_off[i].var
		<< ";" << setprecision(5) << stat_ns_thr[i].var
		<< ";" << setprecision(5) << stat_ns_sen[i].var

		<< ";" << setprecision(5) << jstat_start_offsen[i].cov 
		<< ";" << setprecision(5) << (jstat_start_offsen[i].cov == 0 ? 0 : jstat_start_offsen[i].cov / (sqrt(stat_start_off[i].var) * sqrt(stat_start_sen[i].var)))

		<< ";" << setprecision(5) << jstat_start_offthr[i].cov 
		<< ";" << setprecision(5) << (jstat_start_offthr[i].cov == 0 ? 0 : jstat_start_offthr[i].cov / (sqrt(stat_start_off[i].var) * sqrt(stat_start_thr[i].var)))
		
        << ";" << setprecision(5) << jstat_start_thrsen[i].cov 
		<< ";" << setprecision(5) << (jstat_start_thrsen[i].cov == 0 ? 0 : jstat_start_thrsen[i].cov / (sqrt(stat_start_sen[i].var) * sqrt(stat_start_thr[i].var)));
    }

        DataFile << ";" << msurvivors
        << ";" << fsurvivors
		<< ";" << sum_sexes
        << ";" << setprecision(5) << m_matings 
        << ";" << setprecision(5) <<  var_m_matings
        << ";" << setprecision(5) <<  f_matings
        << ";" << setprecision(5) <<  var_f_matings
        << ";" << setprecision(5) << unmated_f 
        << ";" << setprecision(5) << meanfrs 
        << ";" << setprecision(5) << meanmrs 
        << ";" << setprecision(5) << varfrs 
        << ";" << setprecision(5) << varmrs 
        << ";" << endl;
}

void WriteDataHeaders()
{
	DataFile << "generation";

    for (int i = 1; i <= nTraits; ++i)
    {
		DataFile << ";mean_off" << i 
		<< ";mean_thr" << i
		<< ";mean_sen" << i

        << ";ci_off"<< i
        << ";ci_thr"<< i
        << ";ci_sen"<< i

		<< ";mean_phen_off" << i
		<< ";mean_phen_thr" << i
		<< ";mean_phen_sen" << i

		<< ";mean_ns_off" << i
		<< ";mean_ns_thr" << i
		<< ";mean_ns_sen"<< i
		
        << ";mean_ns_off_ci" << i
		<< ";mean_ns_thr_ci" << i
		<< ";mean_ns_sen_ci"<< i

		<< ";var_off" << i
		<< ";var_thr" << i
		<< ";var_sen" << i

		<< ";var_phen_off" << i
		<< ";var_phen_thr" << i
		<< ";var_phen_sen" << i

		<< ";var_ns_off" << i
		<< ";var_ns_thr" << i
		<< ";var_ns_sen" << i

		<< ";cov_offsen" << i
		<< ";corr_offsen" << i
		
        << ";cov_offthr" << i
		<< ";corr_offthr" << i

        << ";cov_thrsen" << i
		<< ";corr_thrsen" << i;
    }

    DataFile << ";surviving_males"
        << ";surviving_females"
		<< ";N"
        << ";mean_matings_per_male"
        << ";var_matings_per_male"
        << ";mean_matings_per_female"
        << ";var_matings_per_female"
        << ";unmated_f"
        << ";meanfrs"
        << ";meanmrs"
        << ";varfrs"
        << ";varmrs;"
		<< endl;
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
