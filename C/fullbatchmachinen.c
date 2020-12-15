#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

// This file finds the optimal batch size for the binary batch macine that
// extracts all of the work in a batch so you can make a
// contour plot. There is a large recion in the middle wher the optimum is
// always 1 so instead of wasting a lot of time doing a detailed sampling of
// that I sample the edges more finely.
int main(int argc, char *argv[])
{
	int asteps = atoi(argv[1]);
	int ksteps = atoi(argv[2]);
	// Remember that because we have to use the previous batch the function
	// calculates the probabilities for 2*maxn
	int maxn = atoi(argv[3]);
	
	double da = 1.0 / asteps;
	// double dk = 1.0 / ksteps;
	double dk = 0.25 / ksteps;
	
	struct opt optimum;
	
	for(double a = 0;a <= 0.5;a += da)
	{
		// This prints to a different stream to standard output
		// This is not redirected to a file but you can see in in the
		// console window
		fprintf(stderr, "a=%f\n",a);
		
		for(double k = dk;k < 1;k += dk)
		{
			// This prints to a different stream to standard output
			// This is not redirected to a file but you can see in in the
			// console window
			fprintf(stderr, "\tk=%f\n",k);
			
			optimum = extractallfrombatchmaxworkpern(a,k,k,maxn);
			
			printf("%f %f %d\n", a, k, optimum.n);
			// a and 1-a are equilvalent so
			printf("%f %f %d\n", 1 - a, k, optimum.n);
		}
		
		
	}


}