#ifndef _R_PERT_H_
#define _R_PERT_H_

// The hadronic R-ratio in perturbation theory
extern "C" double ru_pert_(double &scms, double &mu_msbar, double &alpha_s, double &mass_c, double &mass_b, double &mass_t);
extern "C" double rd_pert_(double &scms, double &mu_msbar, double &alpha_s, double &mass_c, double &mass_b, double &mass_t);
extern "C" double rs_pert_(double &scms, double &mu_msbar, double &alpha_s, double &mass_c, double &mass_b, double &mass_t);
extern "C" double rc_pert_(double &scms, double &mu_msbar, double &alpha_s, double &mass_c, double &mass_b, double &mass_t);
/*
	scms			center of mass energy squared
  mu_msbar  MSbar renormalization scale
	alpha_s		strong coupling constant at the Z boson mass scale
	mass_c		charm quark mass
	mass_b		bottom quark mass
	mass_t		top quark mass

	See parameters.f for further control over the calculation:
	- The choice of scheme (MSbar vs pole), which determines the numerical values of the charm, bottom and top
	 	quark masses that are appropriate to use.
	- Heavy quark decoupling scales for the running of alpha_s
	- Order of the alpha_s expansion and mass expansion
*/

#endif
