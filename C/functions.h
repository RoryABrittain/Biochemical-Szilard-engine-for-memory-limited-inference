#ifndef FUNCTIONS_HEADER
#define FUNCTIONS_HEADER

double pyn(double a, double k0, double k1, int n, int *y);

double entropy(double a, double k0, double k1, int n);

struct opt
{
	int n;
	double optimum;
	int m;
};

void pi(double a, double k0, double k1, int n, double *p);

struct meaninf
{
	double mean;
	double infimum;
};

void pijoint(double a, double k0, double k1, int n, double p[64][64]);

struct opt binarybatchworkusememory(double a, double k0, double k1, int n);

struct opt binarybatchusememorymaxworkpern(double a, double k0, double k1, int maxn);

struct meaninf binarybatchmeaninfimumusememory(double a, double k0, double k1, int molecules, int trajectories);


struct opt simulatebatchworkoptimum(double a, double k0, double k1, int maxn, int len);

double doublefactorial(int n);

double extractallfrombatchwork(double a, double k0, double k1, int n);

struct opt extractallfrombatchmaxworkpern(double a, double k0, double k1, int maxn);

#endif