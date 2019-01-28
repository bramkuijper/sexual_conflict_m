// RANDOM.CPP
// Truly random number generator
// Obtained from Sander van Doorn
// Added this text at 12-06-2003

#include <math.h>
#include <cstdlib>
#include <iostream>
#include <limits.h>

long int idum;
long int idum2;
#define MBIG LONG_MAX
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
double Uniform(void)
  {
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  if (idum < 0 || iff == 0)
    {
    iff=1;
    mj=MSEED-(idum < 0 ? -idum : idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++)
      {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
      }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++)
        {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
        }
    inext=0;
    inextp=31;
    idum=1;
    }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj * FAC;
  } //End of: double Uniform(void)

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

double gasdev(void)
  {
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;
  double Uniform();
  if  (iset == 0)
    {
    do
      {
      v1=2.0*Uniform()-1.0;
      v2=2.0*Uniform()-1.0;
      r=v1*v1+v2*v2;
      }
    while (r >= 1.0);

    fac=sqrt(-2.0*log(r)/r);
    gset=v1*fac;
    iset=1;
    return v2*fac;
    } //End of: if (iset == 0)
    else
    {
    iset=0;
    return gset;
    } //End of: if (iset == 0) ELSE-part
  } //End of: double gasdev(void)

void SetSeed(int seed)
  {
  double x;
  int i;
  idum = long(-abs(seed));
  idum2 = long(abs(seed));
  for (i=0; i<100; ++i) x= Uniform();
  return;
  } //End of: void SetSeed(int seed)


double Normal(double mean,double stdev)
  {
  return gasdev()*stdev + mean;
  } //End of: double Normal(double mean,double stdev)

int RandomNumber(int N)
  {
  double x;
  int out;

  double k = Uniform();

  x=floor(k*N);
  out=int (x);
  return out;
  }  //End of: int RandomNumber(int N)

int Poison(int mean)
{
    double ran=Uniform();
    double P=exp(-mean);
    double sum=P;
    int i;

    if(ran<sum) {return 0;}
    for(i=1;i<mean*5;++i)
    {
        P=P*mean/i;
        sum=sum+P;
        if(ran<sum) {return i;}
    }
    return mean*5;
}

// random function to do + or - 1

#define IB1 1 //Powers of 2.
#define IB2 2
#define IB5 16
#define IB18 131072
#define MASK (IB1+IB2+IB5)

//Returns as an integer a random bit, based on the 18 low-significance bits in iseed (which is
//modified for the next call).
int rand_pos_neg()
{
	if (idum2 & IB18) { //Change all masked bits, shift, and put 1 into bit 1.
		idum2=((idum2 ^ MASK) << 1) | IB1;
		return 1;
	} 
	else 
	{ //Shift and put 0 into bit 1.
		idum2 <<= 1;
		return -1;
	}
}		
		

