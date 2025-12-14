#ifndef __KEbin__
#define __KEbin__

#include "r_ratio.hh"

class KEbin {
    public:
    KEbin(double x, double Q2, double eE_, double nE_);
    ~KEbin();

    // center of bin
    double x_center = 0;
    double Q2_center = 0;
    double nE = 0;
    double eE = 0;
    
    // weighted value if we perform correction 
    double weighted_x = 0;
    double weighted_Q2 = 0;
    double weighted_y = 0;
    double weighted_nu = 0;

    double x_collection = 0;
    double Q2_collection = 0;
    double y_collection = 0;
    double nu_collection = 0;
    double a1_collection = 0;

    // stastics of bin
    double n_raw = 0;
    double n_events = 0;
    double n_gen = 0;
    double n_acp = 0;

    // status -- good to proceed for calculation?
    bool is_averaged = false;
    bool use_weight = false;

    // processed info of bin
    double x = 0;
    double q2 = 0;
    double y = 0;
    double nu = 0;

    double R = 0;
    double R_err = 0;

    double ga2 = 0;
    double p1 = 0;
    double p2 = 0;
    double unc_ap = 0;
    double unc_a1 = 0;
    double unc_a2 = 0;

    double f1 = 0;
    double f2 = 0;
    double f1err = 0;
    double f2err = 0;
    double a1 = 0;
    double g1 = 0;

    // cs calculation
    double cs = 0;
    double v_bin = 0;
    double c_acc = 0;
    double c_bin = 0;

    // save functions
    double f2p = 0;
    double f2n = 0;
    double f2_he3 = 0;
    double g1p = 0;
    double g1n = 0;

    // functions
    void set_bin(double x_, double Q2_,  double n_, double n_raw_);
    void set_bin_weighted_value(double x, double Q2, double y, double nu);
    void set_volume(double v_);
    void fill_xq2(double x, double Q2);
    void fill_bin(double data_x, double data_Q2, double data_y, double data_nu);

    void save_f1(double f, double ferr);
    void save_f2(double f, double ferr);
    void save_functions(double f2p_, double f2n_, double f2_he3_);

    void average_bin();
    double average_collection(double col);
    void set_to_weight();
    
    void set_n_raw(double n_raw_);
    void calc_cs(double lumi);
    void calc_accept();
    void process_bin(double m_nucleon, double epol, double ipol);
};

KEbin::KEbin(double x, double Q2, double eE_, double nE_)
: x_center(x), Q2_center(Q2), eE(eE_), nE(nE_)
{}

KEbin::~KEbin()
{}

void KEbin::set_bin(double x_, double Q2_, double n_, double n_raw_)
{
  x = x_;
  q2 = Q2_;
  n_events = n_;
  // n_raw = n_raw_;
  is_averaged = true;

  return;
}

void KEbin::set_n_raw(double n)
{
  n_raw = n;

  return;
}


void KEbin::set_volume(double v_)
{
  v_bin = v_;

  return;
}

void KEbin::set_bin_weighted_value(double x_, double Q2_, double y_, double nu_)
{
  weighted_x = x_;
  weighted_Q2 = Q2_;
  weighted_y = y_;
  weighted_nu = nu_;
  is_averaged = true;

  return;
}

void KEbin::fill_xq2(double x, double Q2)
{
  x_collection += x;
  Q2_collection += Q2;

  return;
}

void KEbin::fill_bin(double data_x, double data_Q2, double data_y, double data_nu)
{
  x_collection += data_x;
  Q2_collection += data_Q2;
  y_collection += data_y;
  nu_collection += data_nu;

  n_events ++;

  return;
}

double KEbin::average_collection(double col)
{
  if (n_events == 0)
    return 0;

  return col / n_events;
}

void KEbin::average_bin()
{
  if ( n_events == 0 )
    return;

  weighted_x = average_collection(x_collection);
  weighted_Q2 = average_collection(Q2_collection);
  weighted_y = average_collection(y_collection);
  weighted_nu = average_collection(nu_collection);
  // weighted_a1 = average_collection(a1_collection);

  is_averaged = true;

  return;
}

void KEbin::set_to_weight()
{
  x = weighted_x;
  q2 = weighted_Q2;

  return;
}

void KEbin::save_functions(double f2p_, double f2n_, double f2_he3_)
{
  f2p = f2p_;
  f2n = f2n_;
  f2_he3_ = f2_he3;
}

void KEbin::save_f1(double f, double ferr)
{
  f1 = f;
  f1err = ferr;
  return;
}

void KEbin::save_f2(double f, double ferr)
{
  f2 = f;
  f2err = ferr;
  return;
}

void KEbin::calc_accept()
{
  // acceptance
  if ( n_gen != 0 )
    c_acc = n_acp / n_gen;

  if ( n_acp != 0 )
    c_bin = n_raw / n_acp;

  // if ( n_gen != 0 )
    // printf("c_acc %f, c_bin %f, c_tot %f\n", c_acc, c_bin, c_acc*c_bin);

  return;
}

void KEbin::calc_cs(double lumi)
{
  if ( lumi == 0 || n_events == 0)
    return;

  // prepare constants
  double mp = 938.272088*1e-3;
  double me = 0.510998950*1e-3;
  double ep = 275;
  double ee = 18;
  double Q2max = 4*ep*ee;
  double nu_max = Q2max/(2*mp);
  nu = q2/(2*mp*x);
  y = nu / nu_max;
 
  double alpha = 1./137.03603;
  double twoPiAlphaSq = TMath::TwoPi() * alpha * alpha; 
  double hbarcSq = (0.1973269804*0.1973269804)*1.e7; // units of GeV^2 nb

  // double alpha = 0.0072973525643;

  double y_plus = 1 + (1-y)*(1-y);
  // double factor = weighted_Q2 * weighted_Q2 * weighted_x / (2 * TMath::Pi() * alpha * alpha * y_plus);
  // double factor = weighted_Q2 * weighted_Q2 * weighted_x / y_plus;
  double dcs = n_events / (c_acc * c_bin * lumi * v_bin);
  cs = dcs * (q2*q2*x) / (twoPiAlphaSq * y_plus * hbarcSq * 1e6);

  // cs = n_events / lumi * factor;

  return;
}

void KEbin::process_bin(double m_nucleon, double ePol, double iPol)
{
  // if ( !is_averaged )
  //   return;

  if ( use_weight )
    set_to_weight();

  find_R1998(x, q2, R, R_err);
  
  double nM = m_nucleon;
  double eM = 0.000511;
  double nPz = sqrt(nE*nE - nM*nM);
  double ePz = -1*sqrt(eE*eE - eM*eM);
  
  ROOT::Math::PxPyPzEVector l(0, 0, ePz, eE);
  ROOT::Math::PxPyPzEVector p(0, 0, nPz, nE);

  // TLorentzVector fix_p(0, 0, nPz, nE);
  // TLorentzVector fix_l(0, 0, ePz, eE);
  // TLorentzVector boost(-fix_p[0],-fix_p[1],-fix_p[2],fix_p[3]);
  // TVector3 b = boost.BoostVector();
  // fix_p.Boost(b);
  // fix_l.Boost(b);
  // cout << Form("(%f, %f, %f, %f)", fix_l[0], fix_l[1], fix_l[2], fix_l[3]) << Form(" (%f, %f, %f, %f)", fix_l[0], fix_l[1], fix_l[2], fix_l[3]) << endl;

  double s = (l + p).M2();
  double nu_max = s / (2*m_nucleon);

  nu = q2/(2*m_nucleon*x);
  y = nu / nu_max;

  /*
  double Eprime = nu_max - nu; // only in fixed target frame
  double s2t = q2 / (4*nu_max*Eprime);
  double t2t = 1. / (1-s2t) - 1;
  double epsilon = 1. / (1+2*(1+nu*nu/q2)*t2t);
  double D = (1 - Eprime*epsilon/nu_max) / (1 + epsilon*R);
  double eta = epsilon*sqrt(q2) / (nu_max - Eprime*epsilon);
  double d = D*sqrt(2*epsilon/(1+epsilon));
  double zeta = eta*(1+epsilon) / (2*epsilon);

  double w2  = (m_nucleon*m_nucleon + (2.*m_nucleon*nu) - q2);
  
  // y = nu / eE; // only in fixed target frame
  y = (w2 + q2 - m_nucleon*m_nucleon) / (s -  m_nucleon*m_nucleon);
  cout << eE << " " << nE << " " << x << " " << q2 << " " << nu << " " << y << " " << m_nucleon << endl;
  cout << "var " << s2t << " " << t2t << " " << epsilon << " " << D << " " << eta << " " << d << " " << zeta << endl;

  // std::cout << weighted_y << " " << y << " " << weighted_nu << " " << nu << std::endl;

  double xi = zeta;
  p1 = 1 / (D*(1+eta*xi));
  p2 = eta / (d*(1+eta*xi));

  if ( n_events != 0 )
  {
    unc_ap = 1./(sqrt(n_events/2.) * ePol * iPol);
    unc_a1 = sqrt(pow(p2*unc_ap, 2) + pow(p1*unc_ap, 2));
  }

  cout << "a1 err " << unc_ap << " " << unc_a1 << endl;
  */
  // long double ga2 = 4*m_nucleon*m_nucleon*x*x / q2;
  // long double D = (y*(2-y)*(2+ga2*y)) / (2*(1+ga2)*y*y + (4*(1-y)-ga2*y*y)*(1+R));
  // long double d = sqrt(4*(1-y) - ga2*y*y) / (2-y) * D;
  // long double eta = sqrt(ga2) * (4*(1-y) - ga2*y*y) / ((2-y) * (2+ga2*y));
  // long double xi = sqrt(ga2) * (2-y) / (2+ga2*y);

  ga2 = 4*m_nucleon*m_nucleon*x*x / q2;
  double D = (y*(2-y)*(2+ga2*y)) / (2*(1+ga2)*y*y + (4*(1-y)-ga2*y*y)*(1+R));
  double d = sqrt(4*(1-y) - ga2*y*y) / (2-y) * D;
  double eta = sqrt(ga2) * (4*(1-y) - ga2*y*y) / ((2-y) * (2+ga2*y));
  double xi = sqrt(ga2) * (2-y) / (2+ga2*y);

  p1 = 1 / (D*(1+eta*xi));
  p2 = eta / (d*(1+eta*xi));

  double a2p1 = xi / (D*(1+eta*xi));
  double a2p2 = 1 / (d*(1+eta*xi));

  if ( n_events != 0 )
  {
    unc_ap = 1./(sqrt(n_events/2.) * ePol * iPol);
    unc_a1 = sqrt(pow(p1*unc_ap, 2) + pow(p2*unc_ap, 2));
    unc_a2 = sqrt(pow(a2p1*unc_ap, 2) + pow(a2p2*unc_ap, 2));
  }

  // cout << "a1 err " << unc_ap << " " << unc_a1 << endl;

  // std::cout << "processing .. " << n_events << " " << unc_ap << " " << unc_a1 << std::endl;

  // printf("x %f Q2 %f y %f p1 %f p2 %f unc_ap %f, unc_a1 %f\n", weighted_x, weighted_Q2, y, p1, p2, unc_ap, unc_a1);
}

#endif