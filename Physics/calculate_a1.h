#include <iostream>
#include <cmath>

double find_a1p(double x, double q2)
{
	double alpha = 0.8126;
	double a = 1.2307;
	double b = -0.4128;
	double beta = 0.0303;

	double alpha_err = 0.0488;
	double a_err = 0.1224;
	double b_err = 0.2162;
	double beta_err = 0.1235;

	double a1p = pow(x, alpha)*(a + b*x)*(1 + beta/q2);

	return a1p;
}

double find_a1n(double x, double q2)
{
	double a_n = -0.0490;
	double b_n = -0.1618;
	double c_n = 0.6979;
	double beta_n = 0.7510;

	double a1n = (a_n + b_n*x + c_n*x*x)*(1+ beta_n/q2);

	return a1n;
}