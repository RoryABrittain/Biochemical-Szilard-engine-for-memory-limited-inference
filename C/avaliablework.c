#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"


// This calculates an approximation to the available work for varying a
// and k to make a contour plot
int main(int argc, char *argv[])
{
	int asteps = atoi(argv[1]);
	int ksteps = atoi(argv[2]);
	int maxn = atoi(argv[3]);
	
	double da = 1.0 / asteps;
	double dk = 1.0 / ksteps;
	
	double l2 = log(2);
	
	double avaliablework;
	
	for(double a = 0;a <= 0.5;a += da)
	{
		// This prints to a different stream to standard output
		// This is not redirected to a file but you can see in in the
		// console window
		fprintf(stderr, "a=%f\n",a);

		double k = 0.0000000001;
		avaliablework = l2 - entropy(a, k, k, maxn) / maxn;
		
		printf("%f %f %.16f\n", a, k, avaliablework);
			// a and 1-a are equilvalent so
		printf("%f %f %.16f\n", 1 - a, k, avaliablework);
		
		k = 1;
		avaliablework = l2 - entropy(a, k, k, maxn) / maxn;
		
		printf("%f %f %.16f\n", a, k, avaliablework);
			// a and 1-a are equilvalent so
		printf("%f %f %.16f\n", 1 - a, k, avaliablework);
		
		for(double k = dk;k < 1.0;k += dk)
		{
			avaliablework = l2 - entropy(a, k, k, maxn) / maxn;
			
			printf("%f %f %.16f\n", a, k, avaliablework);
			// a and 1-a are equilvalent so
			printf("%f %f %.16f\n", 1 - a, k, avaliablework);
		}
	}
	
}