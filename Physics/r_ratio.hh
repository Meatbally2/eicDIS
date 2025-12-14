#include <math.h>

void find_R1998(double x, double Q2, double &R, double &dR)
{
	// https://arxiv.org/pdf/hep-ex/9808028.pdf

	double theta = 1 + 12 * (Q2/(Q2+1)) * ((0.125*0.125)/(0.125*0.125+x*x));

	// 

	double a1 = 0.0485;
	double a2 = 0.5470;
	double a3 = 2.0621;
	double a4 = -0.3804;
	double a5 = 0.5090;
	double a6 = -0.0285;

	double Ra = a1 / log(Q2/0.04) * theta + (a2 / pow(pow(Q2,4)+pow(a3,4), 1/4.)) * (1 + a4*x + a5*x*x) * pow(x, a6);

	//

	double b1 = 0.0481;
	double b2 = 0.6114;
	double b3 = -0.3509;
	double b4 = -0.4611;
	double b5 = 0.7172;
	double b6 = -0.0317;

	double Rb = b1 / log(Q2/0.04) * theta + ((b2/Q2) + (b3/(Q2*Q2+0.3*0.3))) * (1 + b4*x +b5*x*x) * pow(x, b6);

	//

	double c1 = 0.0577;
	double c2 = 0.4644;
	double c3 = 1.8288;
	double c4 = 12.3708;
	double c5 = -43.1043;
	double c6 = 41.7415;

	double Q2thr = c4*x + c5*x*x + c6*x*x*x;

	double Rc = c1 / log(Q2/0.04) * theta + c2 * pow((Q2 - Q2thr)*(Q2 - Q2thr) + c3*c3, -1/2.);

	// Average of three models

	R = (Ra + Rb + Rc) / 3.;

	// 

	dR = 0.0078 - 0.013*x + (0.070 - 0.39*x + 0.70*x*x) / (1.7 + Q2);

	return;
}


void find_R1990(double x, double Q2, double &R, double &dR)
{
	// https://www.slac.stanford.edu/pubs/slacreports/reports11/slac-r-357.pdf
	
	double theta = 1 + 12 * (Q2/(Q2+1)) * ((0.125*0.125)/(0.125*0.125+x*x));

	// 

	double a1 = 0.0672;
	double a2 = 0.4671;
	double a3 = 1.8979;

	double Ra = a1 / log(Q2/0.04) * theta + a2 / pow(pow(Q2,4)+pow(a3,4), 1/4.);

	//

	double b1 = 0.0635;
	double b2 = 0.5747;
	double b3 = -0.3534;

	double Rb = b1 / log(Q2/0.04) * theta + (b2/Q2) +(b3/(Q2*Q2+0.3*0.3));

	//

	double c1 = 0.0599;
	double c2 = 0.5088;
	double c3 = 2.1081;

	double Q2thr = 5 * pow(1-x, 5);

	double Rc = c1 / log(Q2/0.04) * theta + c2 * pow((Q2 - Q2thr)*(Q2 - Q2thr) + c3*c3, -1/2.);

	// Average of three models

	R = (Ra + Rb + Rc) / 3.;

	// 

	double B = fmax(0.05, 8.33 * x - 0.66);
	double dRlow = 0.020 + (0.006 + 0.03 * x * x) * abs(log(Q2/B));
	double dRhigh = (0.1 * pow(x,20)) / (pow(0.86,20) + pow(x,20));
	double dRmodle = pow((Ra-R)*(Ra-R)/2. + (Rb-R)*(Rb-R)/2. + (Rc-R)*(Rc-R)/2., 1/2.);

	dR = pow(dRlow*dRlow + dRhigh*dRhigh + dRmodle*dRmodle, 1/2.);

	return;
}

void calc_r_ratio(double x, double Q2, double &r_ratio, double &r_error)
{
	if ( x < 0.03 )
	{
		// fit to 0.005 < x < 0.86, 0.5 < Q2 < 130 (GeV/c)^2 
		find_R1998(x, Q2, r_ratio, r_error);
	}
	else
	{
		// fit to 0.1 < x < 0.9, 0.6 < Q2 < 20 (GeV/c)^2 
		// not supposed to be use for Q2 < 0.3 (GeV/c)^2
		find_R1990(x, Q2, r_ratio, r_error);
	}

	return;
}