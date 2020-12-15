#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "functions.h"


double pyn(double a, double k0, double k1, int n, int *y)
{

  double vector[2] = {k1/(k0+k1),k0/(k0+k1)};

  int i;
  double temp[2];

  for(i=0;i<n;i++)
  {
    if(*(y+i))
    {
      temp[0] = a*((1-k0)*vector[0]+k1*vector[1]);
      temp[1] = (1-a)*(k0*vector[0]+(1-k1)*vector[1]);
    }
    else
    {
      temp[0] = (1-a)*((1-k0)*vector[0]+k1*vector[1]);
      temp[1] = a*(k0*vector[0]+(1-k1)*vector[1]);
    }

    vector[0] = temp[0];
    vector[1] = temp[1];
  }

  double probability = vector[0]+vector[1];
  return(probability);
}

double entropy(double a, double k0, double k1, int n)
{

	//unsigned long long goes up to 2^64-1 so that limits n
	//because 1<<64=0. I need an extra condition so it loops through
	//all possible values of an unsigned long long in that case
	if(n>63)
	{
		printf("Error: Cannot calculate entropy for n>63\n");
		return(NAN);
	}
	int y[64];

	double probability;
	double entropy = 0;
	
	// 1<<n calculates 2^n.
	unsigned long long maxnumber = 1<<n;

	//This is probably not the most efficient way because I don't need to set
	//each value of y each time
	for(unsigned long long i = 0; i < maxnumber; ++i)
	{
		for(int j=0;j<n;j++)
		{
			y[j] = (i>>j)&1;

		}
		
		probability = pyn(a,k0,k1,n,y);
		entropy -= probability * log(probability);
	}

	return(entropy);
}

void pi(double a, double k0, double k1, int n, double *p)
{
	// This function has the side effect of setting the array p equal
	// to the probabilty distribution
	// It assumes that p is initialised to zeros
	// p must have length n+1
	
	//unsigned long long goes up to 2^64-1 so that limits n
	//because 1<<64=0. I need an extra condition so it loops through
	//all possible values of an unsigned long long in that case
	if(n>63)
	{
		printf("Error: Cannot calculate batch work for n>63\n");
		return;
	}
	int y[64];
	
	int countys;
	
	// 1<<n calculates 2^n.
	unsigned long long maxnumber = 1<<n;

	//This is probably not the most efficient way because I don't need to set
	//each value of y each time
	for(unsigned long long i = 0; i < maxnumber; ++i)
	{
		countys = 0;
		
		for(int j=0;j<n;j++)
		{
			y[j] = (i>>j)&1;
			countys += y[j];
		}
		p[countys] += pyn(a,k0,k1,n,y);
	}

	return;
}

void pijoint(double a, double k0, double k1, int n, double p[64][64])
{
	// This function has the side effect of setting the array p equal
	// to the joint probability distribution of  the number of Y* molecules 
	// in 2 successive batches
	// It assumes that p is initialised to zeros
	// p must have length (n+1)*(n+1)
	
	//unsigned long long goes up to 2^64-1 so that limits n
	//because 1<<64=0. I need an extra condition so it loops through
	//all possible values of an unsigned long long in that case
	if(2 * n>63)
	{
		printf("Error: Cannot calculate batch work for n>63\n");
		return;
	}
	int y[64];
	
	int countys1;
	int countys2;
	
	// 1<<n calculates 2^n.
	unsigned long long maxnumber = (unsigned long long)1<<(2 * n);

	//This is probably not the most efficient way because I don't need to set
	//each value of y each time
	for(unsigned long long i = 0; i < maxnumber; ++i)
	{
		countys1 = 0;
		countys2 = 0;
		
		for(int j=0;j<n;j++)
		{
			y[j] = (i>>j)&1;
			countys1 += y[j];
		}
		
		for(int j=n;j<(2 * n);j++)
		{
			y[j] = (i>>j)&1;
			countys2 += y[j];
		}
		
		// for(int i = 0;i < 2 * n;i++)
		// {
			// printf("%d",y[i]);
		// }
		// printf("\n");
		
		p[countys1][countys2] += pyn(a,k0,k1,2 * n,y);
	}
	
	// for(int i = 0;i < n;i++)
	// {
		// for(int j = 0;j < n;j++)
		// {
			// if(p[i][j] > 0)
			// {
				// printf("%d %d %f\n", i, j, p[i][j]);
			// }
		// }
	// }

	return;
}

struct opt binarybatchworkusememory(double a, double k0, double k1, int n)
{
	// This function returns the work extracted by the binary batch machine
	// when the machine uses its previous measurement to make the new
	// measurement
	// This optimises over m (the place where you split the possible
	// batches)
	
	// This is the probability distribution for the number of Y* molecules
	// in a batch
	double p[64];
	
	for(int i=0;i<64;i++)
	{
		p[i] = 0;
	}
	
	pi(a, k0, k1, n, p);
	
	// We need to find the joint distribution for two successive batches
	// we can find the correlation between successive measurements
	double pjoint[64][64];
	
	for(int i=0;i<64;i++)
	{
		for(int j=0;j<64;j++)
		{
			pjoint[i][j] = 0;
		}
	}
	
	pijoint(a, k0, k1, n, pjoint);
	
	struct opt max;
	// If m=0 or N+1 then the work is the single site work so it can't be
	// less than zero
	max.optimum = 0;
	
	// If N=1 the later code doesn't work but I can just return the single
	// site work here
	// in this case 'single site' means the protocol optimal for a Markov
	// input
	if(n==1)
	{
		// log(2) - H[Y2|Y1]
		max.optimum = log(2) + pjoint[0][0] * log(pjoint[0][0] / p[0])
		                     + pjoint[0][1] * log(pjoint[0][1] / p[0])
							 + pjoint[1][0] * log(pjoint[1][0] / p[1])
							 + pjoint[1][1] * log(pjoint[1][1] / p[1]);
		max.m = 1;
		return(max);
	}
	
	// This is 65 because you need to be able to have m=N+1 so you have
	// all of the possilbe states in the first set.
	// Of course, I won't be calculating n=64 anyway
	double work[65];
	
	double l2 = log(2);
	
	// p_in(i<m)
	double pin1;
	
	// p_in(i>=m)
	double pin2;
	
	// expected i/N given i<m
	double p1;
	
	// expected i/N given i>=m
	double p2;
	
	// p(X_{i+1}<m,X_i<m)
	double pll;
	
	// p(X_{i+1}<m,X_i>=m)
	double plm;
	
	// p(X_{i+1}>=m,X_i<m)
	double pml;
	
	// p(X_{i+1}<m,X_i<m)
	double pmm;
	
	// Because of the problems caused by 0log0 I need to consider some m
	// as special cases
	
	// m=0 or n+1 gives the single site machine work
	// In this case you don't do any measurement so I haven't changed it
	// for the case when the machine uses the previous measurement
	// the N=1 case should be better than this so it shouldn't matter anyway
	p1 = 0;
		
	for(int i=0;i<=n;i++)
	{
		p1 += p[i] * i / n;
	}
	
	work[0] = work[n+1] = n * (l2 + p1 * log(p1) + (1 - p1) * log(1 - p1));
	
	// m=1 means p1=0
	pin1 = p[0];
		
	pin2 = 0;
		
	for(int i=1;i<=+n;i++)
	{
		pin2 += p[i];
	}
		
	p2 = 0;
		
	for(int i=1;i<=n;i++)
	{
		p2 += p[i] / pin2 * i / n;
	}
	
	pll = 0;
	for(int i=0;i<1;i++)
	{
		for(int j=0;j<1;j++)
		{
			pll += pjoint[i][j];
		}
	}
	
	plm = 0;
	for(int i=0;i<1;i++)
	{
		for(int j=1;j<=n;j++)
		{
			plm += pjoint[i][j];
		}
	}
	
	pml = 0;
	for(int i=1;i<=n;i++)
	{
		for(int j=0;j<1;j++)
		{
			pml += pjoint[i][j];
		}
	}
	
	pmm = 0;
	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			pmm += pjoint[i][j];
		}
	}
	
	work[1] = n * (l2 + pin2 * (p2 * log(p2) + (1 - p2) * log(1 - p2)))
		       + pll * log(pll / pin1) + plm * log(plm / pin1)
			   + pml * log(pml / pin2) + pmm * log(pmm / pin2);
						
	// m=n means p2=1
	pin1 = 0;
		
	for(int i=0;i<n;i++)
	{
		pin1 += p[i];
	}
		
	pin2 = p[n];
		
	p1 = 0;
		
	for(int i=0;i<n;i++)
	{
		p1 += p[i] / pin1 * i / n;
	}
		
	pll = 0;
	for(int i=0;i<n;i++)
	{
		for(int j=0;j<n;j++)
		{
			pll += pjoint[i][j];
		}
	}
	
	plm = 0;
	for(int i=0;i<n;i++)
	{
		for(int j=n;j<=n;j++)
		{
			plm += pjoint[i][j];
		}
	}
	
	pml = 0;
	for(int i=n;i<=n;i++)
	{
		for(int j=0;j<n;j++)
		{
			pml += pjoint[i][j];
		}
	}
	
	pmm = 0;
	for(int i=n;i<=n;i++)
	{
		for(int j=n;j<=n;j++)
		{
			pmm += pjoint[i][j];
		}
	}
		
	work[n] = n * (l2 + pin1 * (p1 * log(p1) + (1 - p1) * log(1 - p1)))
		       + pll * log(pll / pin1) + plm * log(plm / pin1)
			   + pml * log(pml / pin2) + pmm * log(pmm / pin2);
			   // I think i have actuyally done H[J_i|J_i+1] but it does't matter
			   // which way around it is.
	
	// now I can loop over the other m that don't have a problem
	for(int m=2;m<=n-1;m++)
	{
		pin1 = 0;
		
		for(int i=0;i<m;i++)
		{
			pin1 += p[i];
		}
				
		pin2 = 0;
		
		for(int i=m;i<=+n;i++)
		{
			pin2 += p[i];
		}
				
		p1 = 0;
		
		for(int i=0;i<m;i++)
		{
			p1 += p[i] / pin1 * i / n;
		}
				
		p2 = 0;
		
		for(int i=m;i<=n;i++)
		{
			p2 += p[i] / pin2 * i / n;
		}
		
		pll = 0;
		for(int i=0;i<m;i++)
		{
			for(int j=0;j<m;j++)
			{
				pll += pjoint[i][j];
			}
		}
	
		plm = 0;
		for(int i=0;i<m;i++)
		{
			for(int j=m;j<=n;j++)
			{
				plm += pjoint[i][j];
			}
		}
	
		pml = 0;
		for(int i=m;i<=n;i++)
		{
			for(int j=0;j<m;j++)
			{
				pml += pjoint[i][j];
			}
		}
	
		pmm = 0;
		for(int i=m;i<=n;i++)
		{
			for(int j=m;j<=n;j++)
			{
				pmm += pjoint[i][j];
			}
		}
				
		work[m] = n * (l2 + pin1 * (p1 * log(p1) + (1 - p1) * log(1 - p1))
		                  + pin2 * (p2 * log(p2) + (1 - p2) * log(1 - p2)))
		          + pll * log(pll / pin1) + plm * log(plm / pin1)
			      + pml * log(pml / pin2) + pmm * log(pmm / pin2);
		
		//printf("m=%d, p1=%f, p2=%f, pll=%f, plm=%f, pml=%f, pmm=%f\n", m, p1, p2, pll, plm, pml, pmm);
	}
	
	for(int m=0;m<=n+1;m++)
	{
				  
		if(work[m] > max.optimum)
		{
			max.optimum = work[m];
			max.m = m;
		}
		
	}
	
	return(max);
	
	
}

struct opt binarybatchusememorymaxworkpern(double a, double k0, double k1, int maxn)
{
	// This function maximises the work of the binary batch machine
	// over n when it uses the previous mesurement to make the next
	
	struct opt workpern;
	
	struct opt nandmaxworkpern;
	
	nandmaxworkpern.optimum = 0;
	nandmaxworkpern.n = 1;
	
	for(int n=1;n<=maxn;n++)
	{
		workpern = binarybatchworkusememory(a, k0, k1, n);
		workpern.optimum = workpern.optimum / n;
		if(workpern.optimum > nandmaxworkpern.optimum)
		{
			nandmaxworkpern.optimum = workpern.optimum;
			nandmaxworkpern.m = workpern.m;
			nandmaxworkpern.n = n;
		}
	}
	
	return(nandmaxworkpern);
}

struct meaninf binarybatchmeaninfimumusememory(double a, double k0, double k1, int molecules, int trajectories)
{
	double l2 = log(2);
	
	// Stationary probability that s is 0
	double ps = k1/(k0+k1);
	
	// Set the random seed. This only changes once per second and the
	// first rand() seems to incriment
	unsigned int t = time(NULL);
	srand(t);
	
	// If state is 0 then p(Y^*)=alpha
	int s;
	
	// y=0 is Y^*
	int y;
	
	// p(switch|s_i)
	double sswitch[2] = {k0, k1};
	
	// p(y_i=0|s_i)
	double py[2] = {a, 1-a};
	
	// I need to find the optimum batch size and position to split the
	// batches for the parameters so I know what protocol is being used.
	const int maxbatchsize = 10;
	// It is important to remember to use the correct function here
	struct opt optimal = binarybatchusememorymaxworkpern(a, k0, k1, maxbatchsize);
	// printf("optimal work per molecule=%f\n", optimal.optimum);
	// printf("optimal work in trajectory=%f\n", molecules * optimal.optimum);
	// optimal batch size
	int n = optimal.n;
	// optimal position to split the batches
	int m =optimal.m;
	
	// printf("n=%d, m=%d\n",n,m);
	
	struct meaninf output;

	// if the optimal n is the maximum batch size that I have calculated
	// then the true optimal batch size could be bigger and I wouldn't
	// know
	if(n == maxbatchsize)
	{
		fprintf(stderr,"Greater than maximum batch size\n");
		output.mean = NAN;
		output.infimum = NAN;
		return(output);
	}
	
	// I need the probability distribution to get the resetting work
	double pin[maxbatchsize+1];
	for(int i=0;i<maxbatchsize+1;i++)
	{
		pin[i] = 0;
	}
	pi(a, k0, k1, n, pin);
	
	// Calculate the values will need to calculate the work
	// I don't know but I think it is faster to only take the log once
	
	// p(i<m)
	double pless = 0;
	for(int i=0;i<m;i++)
	{
		pless += pin[i];
	}
	
	// p(i>=m)
	double pmore = 0;
	for(int i=m;i<=n;i++)
	{
		pmore += pin[i];
	}
	
	// p1
	double p1 = 0;
	for(int i=0;i<m;i++)
	{
		p1 += pin[i] / pless * i / n;
		// printf("p1=%f\n",p1);
	}
	double lp1 = log(p1);
	double lmp1 = log(1 - p1);
	
	
	// p2
	double p2 = 0;
	for(int i=m;i<=n;i++)
	{
		p2 += pin[i] / pmore * i / n;
		// printf("p2=%f\n",p2);
	}
	double lp2 = log(p2);
	double lmp2 = log(1 - p2);
	
	// printf("p1=%f\n",p1);
	// printf("p2=%f\n",p2);
	
	// We need to find the joint distribution for two successive batches
	// we can find the correlation between successive measurements
	double pjoint[64][64];
	
	for(int i=0;i<64;i++)
	{
		for(int j=0;j<64;j++)
		{
			pjoint[i][j] = 0;
		}
	}
	
	pijoint(a, k0, k1, n, pjoint);
	
	// this is the condional probability of next measurement based on the
	// previous
	// 0 is less than m and 1 is more than or equal to
	// pgiven[a][b] is a given the previous was b
	// the order is opposite to pjoint which has the first batch first
	double pgiven[2][2] = {0};
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<m;j++)
		{
			pgiven[0][0] += pjoint[i][j];
		}
	}
	pgiven[0][0] /= pless;
	
	for(int i=m;i<=n;i++)
	{
		for(int j=0;j<m;j++)
		{
			pgiven[0][1] += pjoint[i][j];
		}
	}
	pgiven[0][1] /= pmore;
	
	for(int i=0;i<m;i++)
	{
		for(int j=m;j<=n;j++)
		{
			pgiven[1][0] += pjoint[i][j];
		}
	}
	pgiven[1][0] /= pless;
	
	for(int i=m;i<=n;i++)
	{
		for(int j=m;j<=n;j++)
		{
			pgiven[1][1] += pjoint[i][j];
		}
	}
	pgiven[1][1] /= pmore;
	
	// printf("pmore=%f\n",pmore);
	// printf("pless=%f\n",pless);
	
	// ln(p(i|j))
	double lpgiven[2][2];
	for(int i=0;i<2;i++)
	{
		for(int j=0;j<2;j++)
		{
			lpgiven[i][j] = log(pgiven[i][j]);
			// printf("p(%d|%d)=%f, log(p(%d|%d))=%f\n",i,j,pgiven[i][j],i,j,lpgiven[i][j]);
		}
	}
	
	// Initialies things I will need in the loop
	double work;
	
	double infimum;
	
	double averageinfimum = 0;
	
	double averagework = 0;
	
	int countbatch;
	int ys;
	// I also need to store the state of the previous measurement.
	int previousmeasurement;
	
	for(int i=0;i<trajectories;i++)
	{
		work = 0;
		infimum = 0;
		
		// set up first molecule
		s = rand() < ps*RAND_MAX ? 0 : 1;

		
		y = rand() < py[s]*RAND_MAX ? 0 : 1;
		
		ys = 1-y;
		
		// simulate the rest of one batch so I have a previous ys
		for(int j=0;j<n-1;j++)
		{
			if(rand() < sswitch[s]*RAND_MAX)
			{
				// XOR with 1, i.e. switch state
				s = s^1;
			}
			
			y = rand() < py[s]*RAND_MAX ? 0 : 1;
			
			ys += 1-y;
		}
		
		previousmeasurement = ys < m ? 0 : 1;
		
		countbatch = 0;
		ys = 0;
		
		// simulate the actual run of molecules for the machine to extract
		// work from
		for(int j=0;j<molecules;j++)
		{
			if(rand() < sswitch[s]*RAND_MAX)
			{
				// XOR with 1, i.e. switch state
				s = s^1;
			}
			
			
			y = rand() < py[s]*RAND_MAX ? 0 : 1;
			
			countbatch += 1;
			ys += 1-y;
			
			if(countbatch%n==0)
			{
				// When calculating the work remember that what I usually
				// call `i' is called `ys' here.
				
				// I have to do the different cases otherwise it will try
				// to calculate 0log0
				
				if(m==0) // p(i<m) = 0. bias machine only
				{
					work += n * (l2 + (double)ys / n * lp2
					                + (1 - (double)ys / n) * lmp2);
				}
				else if(n==1 && m==1)
				{
					// printf("prev=%d, ys=%d, p(i|j)=%f\n",previousmeasurement,ys,pgiven[ys][previousmeasurement]);
					// printf("log(p(i|j))=%f\n",lpgiven[ys][previousmeasurement]);
					work += l2 + lpgiven[ys][previousmeasurement];
				}
				else if(m==1) // p1 = 0
				{
					if(ys < m)
					{
						work += n * l2 + lpgiven[0][previousmeasurement];
					}
					else
					{
						work += n * (l2 + (double)ys / n * lp2
					                    + (1 - (double)ys / n) * lmp2)
								+ lpgiven[1][previousmeasurement];
					}
				}		
				else if(m==n) // p2 = 1
				{
					if(ys < m)
					{
						work += n * (l2 + (double)ys / n * lp1
					                    + (1 - (double)ys / n) * lmp1)
								+ lpgiven[0][previousmeasurement];
					}
					else
					{
						work += n * l2 + lpgiven[1][previousmeasurement];
					}
				}
					
				else if(m==n + 1) // p(i>=m) = 0
				{
					work += n * (l2 + (double)ys / n * lp1
					                + (1 - (double)ys / n) * lmp1);
				}
				else
				{
					if(ys < m)
					{
						work += n * (l2 + (double)ys / n * lp1
				                        + (1 - (double)ys / n) * lmp1)
								+ lpgiven[0][previousmeasurement];
					}
					else
					{
						work += n * (l2 + (double)ys / n * lp2
				                        + (1 - (double)ys / n) * lmp2)
								+ lpgiven[1][previousmeasurement];
					}
				}
				
				previousmeasurement = ys < m ? 0 : 1;
				ys = 0;
			}
			
		
			if(work < infimum)
			{
				infimum = work;
			}
			
			// printf("trajectory=%d, molecule=%d, work=%f\n",i,j,work);
			
		}
		
		averageinfimum += infimum;
		averagework += work;
		
	}
	
	averageinfimum /= trajectories;
	averagework /= trajectories;
	
	output.mean = averagework;
	output.infimum = averageinfimum;
	
	return(output);
}

struct opt simulatebatchworkoptimum(double a, double k0, double k1, int maxn, int len)
{
	// This funciton might be faster if i don't shift each element of the array
	// that stores the state of he moleculse allong on in each step.
	
	// we need the maximum possible batch size to create the variables to
	// store the results with enough space
	enum maxbatchsize {maxbatchsize = 1000};
	
	// in this version we only remember the last n inputs
	int data[maxbatchsize];
	
	// Stationary probability that s is 0
	double ps = (k0+k1 == 0) ? 0.5 :k1/(k0+k1);

	
	// Set the random seed. This only changes once per second and the
	// first rand() seems to incriment
	srand(time(NULL));
	
	// If state is 0 then p(Y^*)=alpha
	int s;
	
	// y=0 is Y^*
	int y;
	
	// p(switch|s_i)
	double sswitch[2] = {k0, k1};
	
	// p(y_i=0|s_i)
	double py[2] = {a, 1-a};
	
	// number of Y* in a batch. the element shows the batch size
	int batch[maxbatchsize];
	int prevbatch[maxbatchsize];
	
	// p_in(i<m)
	int pin1count[maxbatchsize] = {0};
	double pin1;
	
	// p_in(i>=m)
	double pin2;
	
	// expected i/N given i<m
	int p1count[maxbatchsize] = {0};
	double p1;
	
	// expected i/N given i>=m
	int p2count[maxbatchsize] = {0};
	double p2;
	
	// p(X_{i+1}<m,X_i<m)
	int pllcount[maxbatchsize] = {0};
	double pll;
	
	// p(X_{i+1}<m,X_i>=m)
	double plm;
	
	// p(X_{i+1}>=m,X_i<m)
	int pmlcount[maxbatchsize] = {0};
	double pml;
	
	// p(X_{i+1}<m,X_i<m)
	double pmm;
	
	double work;
	double l2 = log(2);
	
	// Initialise states
	
	s = rand() < ps*RAND_MAX ? 0 : 1;

		
	y = rand() < py[s]*RAND_MAX ? 0 : 1;
	
	data[0] = y;
	
	// burn in so initially there is a history in data[]
	// this is not strictly necessary but i think it makes the main loop cleaner
	// if you don't have to keep checking that you have had at least one batch
	for(int i = 1;i <= maxbatchsize;i++)
	{
		if(rand() < sswitch[s]*RAND_MAX)
		{
			// XOR with 1, i.e. switch state
			s = s^1;
		}
			
			
		y = rand() < py[s]*RAND_MAX ? 0 : 1;
		
		for(int j = maxbatchsize - 1;j >= 0;j--)
		{
			data[j + 1] = data[j];
		}
		data[0] = y;
		
		for(int n = 1;n <= maxn;n++)
		{
			// check that we have reached the end of a batch
			if(i % n == 0)
			{
				prevbatch[n] = batch[n];
				
				batch[n] = 0;
				for(int j = 0;j < n;j++)
				{
					batch[n] += 1 - data[j]; //because y=0 is Y*
				}
			}
		}
		
	}
	
	// This is the loop over the successive molecule in one trajectory
	// starting at 1 so you don't think you've reached the end of a batch of all
	// sizes on the first molecule
	for(int i = 1;i <= len;i++)
	{
		if(rand() < sswitch[s]*RAND_MAX)
		{
			// XOR with 1, i.e. switch state
			s = s^1;
		}
			
			
		y = rand() < py[s]*RAND_MAX ? 0 : 1;
		// printf("%d", 1 - y);
		
		for(int j = maxbatchsize - 1;j >= 0;j--)
		{
			data[j+1] = data[j];
		}
		data[0] = y;
		
		// Loop over batch sizes
		for(int n = 1;n <= maxn;n++)
		{
			// check that we have reached the end of a batch
			if(i % n == 0)
			{
				// if(n == 5)
				// {
					// printf("\n");
					// for(int l=0;l<5;l++)
					// {
						// printf("%d", 1 - data[l]);
					// }
				// }
				
				prevbatch[n] = batch[n];
				
				batch[n] = 0;
				for(int j = 0;j < n;j++)
				{
					batch[n] += 1 - data[j]; //because y=0 is Y*
					// if(n == 5)
					// {
						// printf("\nbatch=%d", batch[n]);
					// }
				}
				
				// if(n == 5)
				// {
					// printf("\nprevbatch=%d, batch=%d\n", prevbatch[n], batch[n]);
				// }
				
				if(batch[n] <= n / 2) // remember this is integer division
				{
					pin1count[n] += 1;
					p1count[n] += batch[n];
				
					if(prevbatch[n] <= n / 2)
					{
						pllcount[n] += 1;
					}

				}
				else
				{
					p2count[n] += batch[n];
				
					if(prevbatch[n] <= n / 2)
					{
						pmlcount[n] += 1;
					}
				}
				
				
				
				
			}
		}
		
	}
	
	struct opt optimum;
	optimum.optimum = 0;
	optimum.m = -1;
	
	
	// Loop over batch sizes
	for(int n = 1;n <= maxn;n++)
	{
		pin1 = (double)pin1count[n] / (len / n);
		pin2 = 1 - pin1;
		p1 = (double)p1count[n] / (n * pin1count[n]);
		p2 = (double)p2count[n] / (n * ((len / n) - pin1count[n]));
		pll = (double)pllcount[n] / (len / n);
		plm = pin1 - pll;
		pml = (double)pmlcount[n] / (len / n);
		pmm = pin2 - pml;
		
		
		// if(n == 5)
		// {
			// printf(
			// "pin1=%f, pin2=%f, p1=%f, p2=%f, pll=%f, plm=%f, pml=%f, pmm=%f\n",
			// pin1, pin2, p1, p2, pll, plm, pml, pmm);
		// }
		
		switch(n)
		{
			case 1: // p1=0 and p2=1
				work = n * l2
				       + pll * log(pll / pin1) + plm * log(plm / pin2)
				       + pml * log(pml / pin1) + pmm * log(pmm / pin2);
				break;
			
			case 2: // p2=1
				work = n * (l2 + pin1 * (p1 * log(p1) + (1 - p1) * log(1 - p1)))
				   + pll * log(pll / pin1) + plm * log(plm / pin2)
				   + pml * log(pml / pin1) + pmm * log(pmm / pin2);
				break;
			
			default :
				work = n * (l2 + pin1 * (p1 * log(p1) + (1 - p1) * log(1 - p1))
		                       + pin2 * (p2 * log(p2) + (1 - p2) * log(1 - p2)))
				   + pll * log(pll / pin1) + plm * log(plm / pin2)
				   + pml * log(pml / pin1) + pmm * log(pmm / pin2);
		}
		
		
		// printf("n=%d, work/n=%f\n", n, work / n);
		if(work / n > optimum.optimum)
		{
			optimum.optimum = work / n;
			optimum.n = n;
		}

	}
	
	return(optimum);
}

double doublefactorial(int n)
{
	
	double temp = 1;
	
	if(n == 0 || n == 1)
	{
		return(1);
	}
	else
	{
		for(int i = n;i > 1;i--)
		{
			temp *= i;
			// printf("%f, ", temp);
		}
	}
	// printf("\n");

	return(temp);
}

double extractallfrombatchwork(double a, double k0, double k1, int n)
{
	// this function finds the work extracted by the batch machine if it fully
	// extracts all work from the batch and the memory records if the batch has
	// <N/2 X* molecules
	
	// we need the equilibrium distribution to calculate the free energy
	double peq[64];
	
	// i use N! repeatedly
	double nfact = doublefactorial(n);
	
	// i use 2^N repeatedly
	double inv2pown = 1.0 / (1<<n);
	
	for(int i = 0;i <= n;i++)
	{
		peq[i] = inv2pown * nfact / (doublefactorial(i) * doublefactorial(n-i));
	}
	
	// This is the probability distribution for the number of Y* molecules
	// in a batch
	double p[64] = {0};
	
	pi(a, k0, k1, n, p);
		
	// the free energy of the batch
	double fin = 0;
	
	//calculate by KL divergence
	for(int i = 0;i <= n;i++)
	{
		fin += p[i] * log(p[i] / peq[i]);
	}
	
	// We need to find the joint distribution for two successive batches
	// we can find the correlation between successive measurements
	double pjoint[64][64] = {0};
	
	pijoint(a, k0, k1, n, pjoint);
	
	// Now find the joint distribution between sucesive 1 bit measurements of if
	// i<N/2
	// first number is fist batch and second number is sectond batch
	double p00 = 0;
	// need this extra logic to make the <N/2 work correctly for integer division
	for(int i = 0;i < (n % 2 == 0 ? n / 2: n / 2 + 1);i++)
	{
		for(int j = 0;j < (n % 2 == 0 ? n / 2: n / 2 + 1);j++)
		{
			p00 += pjoint[i][j];
		}
	}
	
	double p01 = 0;
	for(int i = 0;i < (n % 2 == 0 ? n / 2: n / 2 + 1);i++)
	{
		for(int j = (n % 2 == 0 ? n / 2: n / 2 + 1);j <= n;j++)
		{
			p01 += pjoint[i][j];
		}
	}
	
	double p10 = 0;
	for(int i = (n % 2 == 0 ? n / 2: n / 2 + 1);i <= n;i++)
	{
		for(int j = 0;j < (n % 2 == 0 ? n / 2: n / 2 + 1);j++)
		{
			p10 += pjoint[i][j];
		}
	}
	
	double p11 = 0;
	for(int i = (n % 2 == 0 ? n / 2: n / 2 + 1);i <= n;i++)
	{
		for(int j = (n % 2 == 0 ? n / 2: n / 2 + 1);j <= n;j++)
		{
			p11 += pjoint[i][j];
		}
	}
	
	
	// mutual information between sucessive measurements
	double info;
	
	info = p00 * log(p00 / ((p00 + p01) * (p00 + p10))) +
	       p01 * log(p01 / ((p00 + p01) * (p01 + p11))) +
		   p10 * log(p10 / ((p10 + p11) * (p00 + p10))) +
		   p11 * log(p11 / ((p11 + p10) * (p01 + p11)));
	
	
	double work = 0;
	work = fin + info;
	
	return(work);
	
}

struct opt extractallfrombatchmaxworkpern(double a, double k0, double k1, int maxn)
{
	// This function maximises the work of the binary batch machine
	// over n when it uses the previous mesurement to make the next
	
	double workpern;
	
	struct opt nandmaxworkpern;
	
	nandmaxworkpern.optimum = 0;
	nandmaxworkpern.n = 1;
	
	for(int n=1;n<=maxn;n++)
	{
		workpern = extractallfrombatchwork(a, k0, k1, n) / n;
		
		if(workpern > nandmaxworkpern.optimum)
		{
			nandmaxworkpern.optimum = workpern;
			nandmaxworkpern.n = n;
		}
	}
	
	return(nandmaxworkpern);
}

