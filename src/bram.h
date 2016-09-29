#ifndef DISTREADER_H_
#define DISTREADER_H_

#include <ctime>
#include <string>
#include <iostream>
#include <sstream>
#include <cmath>
//#include <gsl/gsl_cdf.h>

// convert from integer to string
std::string itos(int j)
{
	std::stringstream s;
	s << j;
	return s.str();
}


// filename function 
std::string create_filename(std::string filename)
{
	time_t timep = time(NULL);

	tm *timestruct = localtime(&timep);

	filename.append("_");
	filename.append(itos(timestruct->tm_mday));
	filename.append("_");
	filename.append(itos(timestruct->tm_mon + 1));
	filename.append("_");
	filename.append(itos(timestruct->tm_year + 1900));
	filename.append("_");

	size_t hour = timestruct->tm_hour;

	if (hour < 10)
		filename.append("0");

	filename.append(itos(hour));

	size_t minutes = timestruct->tm_min;

	if (minutes < 10)
		filename.append("0");

	filename.append(itos(minutes));

	size_t seconds = timestruct->tm_sec;

	if (seconds < 10)
		filename.append("0");

	filename.append(itos(seconds));

	filename.append("_");

	// now we are going to insert microsecs
	// to obtain unique filename
	timespec ts;
	if (0!=clock_gettime(CLOCK_REALTIME,&ts))
	{
		throw "error in getting the real time stamp in nanosecs, quitting!";
	}

	filename.append(itos(ts.tv_nsec));

	return(filename);
}

unsigned get_nanoseconds()
{
	timespec ts;

	if (0!=clock_gettime(CLOCK_REALTIME,&ts))
	{
		throw "error in getting the real time stamp in nanosecs, quitting!";
	}

	return(ts.tv_nsec);
}

int linear_search(const int list[], const int maxsize, int value)
{
    int index = 0;
    int position = -1;
    bool found = false;

    while (index < maxsize && !found)
    {
        if (list[index] == value)
        {
            found = true;
            position = index;
        }

        ++index;
    }

    return(position);
}

struct Stats
{
    double mean;
    double mean_ci;
    double sumsquares;
    double var;
    double sumthirds;
    double skew;
    double sumfourths;
    double kurt;
    double sample;
};

struct JointStats
{
    double cov;
    double sum12;
    double sample;
};

void stat_reset(Stats &results)
{
    results.mean = 0;
    results.mean_ci = 0;
    results.sumsquares = 0;
    results.var = 0;
    results.sumthirds = 0;
    results.skew = 0;
    results.sumfourths = 0;
    results.kurt = 0;
    results.sample = 0;
}


void stat_addval(Stats &results, const double traitvalue)
{
    results.mean += traitvalue;   
    results.sumsquares += pow(traitvalue,2);
    results.sumthirds += pow(traitvalue,3);
    results.sumfourths += pow(traitvalue,4);
    ++results.sample;
}

void jstat_reset(JointStats &results)
{
    results.sum12 = 0;
    results.sample = 0;
}

void jstat_addval(JointStats &results, const double trait1, const double trait2)
{
    results.sum12+=trait1*trait2;
    ++results.sample;
}

void jstat_finalize(JointStats &results, const double mean1, const double mean2)
{
    results.cov = results.sum12 / results.sample - mean1 * mean2;
}

void stat_finalize(Stats &results)
{
    if (results.sample == 0)
    {
        results.mean = 0;
        results.var = 0;
        results.mean_ci = 0;
        results.skew = 0;
        results.kurt = 0;
    }
    else
    { 
        results.mean /= results.sample;
        results.var = results.sumsquares == 0 ? 0 : results.sumsquares / results.sample - pow(results.mean,2);

        results.mean_ci = 0; //results.sample > 1 ? gsl_cdf_tdist_Qinv(0.025,results.sample - 1) * results.var / sqrt(results.sample) : 0;
        results.skew = results.sample > 3 ? results.sumthirds / results.sample - 3 * results.mean * results.sumsquares / results.sample + 2 * pow(results.mean,3) : 0;
        results.kurt = results.sample > 4 ? (double) (results.sumfourths / results.sample - 4 * results.sumthirds / results.sample * results.mean + 6 * pow(results.mean,2) * results.sumsquares / results.sample - 3 * pow(results.mean,4)) : 0;

    }
}

#endif
