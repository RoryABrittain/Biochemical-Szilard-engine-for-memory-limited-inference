#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

// This file finds the work for the binary batch macine that
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
	double maxk = atof(argv[4]);
	int entropyn = atoi(argv[5]);

	
	double da = 1.0 / asteps;
	// double dk = 1.0 / ksteps;
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
		
		for(double k = dk;k < maxk;k += dk)
		{
			// This prints to a different stream to standard output
			// This is not redirected to a file but you can see in in the
			// console window
			fprintf(stderr, "\tk=%f\n",k);
			
			optimum = extractallfrombatchmaxworkpern(a,k,k,maxn);
			
			avaliablework = l2 - entropy(a, k, k, entropyn) / entropyn;
			
			efficiency = optimum.optimum / avaliablework;
			
			printf("%f %f %f\n", a, k, efficiency);
			// a and 1-a are equilvalent so
			printf("%f %f %f\n", 1 - a, k, efficiency);
		}
		
	
		
	}
	

}