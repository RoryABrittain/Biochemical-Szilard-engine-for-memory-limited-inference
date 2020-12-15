#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"


// This simulates the optimum work per N of the batch machine for varying a and k
// to make a contour plot of the optimum N
int main(int argc, char *argv[])
{
	int asteps = atoi(argv[1]);
	int ksteps = atoi(argv[2]);
	int maxn = atoi(argv[3]);
	int len = atoi(argv[4]); // number of steps in simulation
	
	double maxk = 0.08;
	
	double da = 1.0 / asteps;
	double dk = maxk / ksteps;
	
	struct opt optimum;
	
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
			
			printf("%f %f %d\n", a, k, optimum.n);
			// a and 1-a are equilvalent so
			printf("%f %f %d\n", 1 - a, k, optimum.n);
		}
	}
	
}