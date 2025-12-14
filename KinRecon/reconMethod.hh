
// electron method
std::vector<float> calc_elec_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
    float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
    float y = 1. - (E/E_ebeam)*TMath::Sin(theta/2)*TMath::Sin(theta/2);
    float x = Q2/(4*E_ebeam*E_pbeam*y);
    return {x, y, Q2};
  }
  
  // jb method
  std::vector<float> calc_jb_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
    float y = sigma_h/(2*E_ebeam);
    float Q2 = pt_had*pt_had / (1-y);
    float x = Q2/(4*E_ebeam*E_pbeam*y);
    return {x, y, Q2};
  }
  
  // double angle method
  std::vector<float> calc_da_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
    float alpha_h = sigma_h/pt_had;
    float alpha_e = TMath::Tan(theta/2);
    float y = alpha_h / (alpha_e + alpha_h);
    float Q2 = 4*E_ebeam*E_ebeam / (alpha_e * (alpha_h + alpha_e));
    float x = Q2/(4*E_ebeam*E_pbeam*y);
    return {x, y, Q2};
  }
  
  // sigma method
  std::vector<float> calc_sig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
    float y = sigma_h/(sigma_h + E*(1 - TMath::Cos(theta))); 
    float Q2 = E*E*TMath::Sin(theta)*TMath::Sin(theta) / (1-y);
    float x = Q2/(4*E_ebeam*E_pbeam*y);
    return {x, y, Q2};
  }
  
  // e-sigma method
  std::vector<float> calc_esig_method(float E, float theta, float pt_had, float sigma_h, float E_ebeam, float E_pbeam) {
    float Q2  = 2.*E_ebeam*E*(1+TMath::Cos(theta));
    float x = calc_sig_method(E,theta,pt_had,sigma_h,E_ebeam,E_pbeam)[0];
    float y = Q2/(4*E_ebeam*E_pbeam*x);
    return {x, y, Q2};
  }