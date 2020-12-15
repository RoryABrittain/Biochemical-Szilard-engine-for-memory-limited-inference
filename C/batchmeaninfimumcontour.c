#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

// This file calculates the mean infiumum for the batch machine at various
// values of a and k to make a contourplot
int main(int argc, char *argv[])
{
	int asteps = atoi(argv[1]);
	int ksteps = atoi(argv[2]);
	int molecules = atoi(argv[3]);
	int trajectories = atoi(argv[4]);
	double maxk = atof(argv[5]);
	
	double da = 1.0 / asteps;
	double dk = maxk / ksteps;
	
	struct meaninf meanandinfimum;
	
	
	for(double a = 0;a <= 0.5;a += da)
	{
		// This prints to a different stream to standard output
		// This is not redirected to a file but you can see in in the
		// console window
		fprintf(stderr, "a=%f\n",a);
		
		for(double k = dk;k <= maxk;k += dk)
		{
			meanandinfimum = binarybatchmeaninfimumusememory(a,k,k,molecules,trajectories);
			
			printf("%f %f %.10f\n", k, a, meanandinfimum.infimum);
			// a and 1-a are equilvalent so
			printf("%f %f %.10f\n", k, 1 - a, meanandinfimum.infimum);
			
		}
		
	}

}