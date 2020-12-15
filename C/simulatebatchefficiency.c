#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"


// This simulates the optimum work per N of the batch machine for varying a
// and k to make a contour plot of the work
int main(int argc, char *argv[])
{
	int asteps = atoi(argv[1]);
	int ksteps = atoi(argv[2]);
	int maxn = atoi(argv[3]);
	int len = atoi(argv[4]); // number of steps in simulation
	int entropyn = atoi(argv[5]);
	
	double maxk = 0.08;
	
	double da = 1.0 / asteps;
	double dk = maxk / ksteps;
	
	struct opt optimum;
	
	double l2 = log(2);
	
	double avaliablework;
	double efficiency;
	
	for(double a = 0;a <= 0.5;a += da)
	{
		// This prints to a different stream to standard output
		// This is not redirected to a file but you can see in in the
		// console window
		fprintf(stderr, "a=%f\n",a);
		
		for(double k = dk;k <= maxk;k += dk)
		{
			// This prints to a different stream to standard output
			// This is not redirected to a file but you can see in in the
			// console window
			// fprintf(stderr, "\tk=%f\n",k);
			
			optimum = simulatebatchworkoptimum(a, k, k, maxn, len);
			
			avaliablework = l2 - entropy(a, k, k, entropyn) / entropyn;
			
			efficiency = optimum.optimum / avaliablework;
			
						
			printf("%f %f %f\n", a, k, efficiency);
			// a and 1-a are equilvalent so
			printf("%f %f %f\n", 1 - a, k, efficiency);
			
		}
	}
	
}