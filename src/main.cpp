#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "core/rhad_pert.h"

int main() {

/******************************************************************************/

  double scms;
  double mu_msbar;
  double alpha_s;
  double mass_c;
  double mass_b;
  double mass_t;
  /*
    scms			center of mass energy squared
    mu_msbar  MSbar renormalization scale
    alpha_s		strong coupling constant at the Z boson mass scale
    mass_c		charm quark mass
    mass_b		bottom quark mass
    mass_t		top quark mass
  */

  alpha_s = 0.1181;
  mass_c = 1.275;
  mass_b = 4.18;
  mass_t = 160;

  std::ofstream outFile;

  const int NUM_POINTS = 1000;
  const double S_MIN = 4.0;
  const double S_MAX = 1000;

  outFile.open("rpert.txt");
  outFile << std::fixed << std::setprecision(6);

  outFile << "----------------------------------------------------------------------------------------------------------------------\n\n"
          << " The contribution to the hadronic R-ratio from u,d,s,c quark flavors. The columns should be read as follows:\n\n"
          << " s(GeV),  R_u,  R_d,  R_s,  R_c, \n\n"
          << " The following inputs were used (masses are in the MSbar scheme): \n\n"
          << " alpha_s(m_z) = " << alpha_s << "\n"
          << " mass_c(mass_c) = " << mass_c << "\n"
          << " mass_b(mass_b) = " << mass_b << "\n"
          << " mass_t(mass_t) = " << mass_t << "\n\n"
          << "----------------------------------------------------------------------------------------------------------------------\n\n";

/******************************************************************************/

  for (int i=0; i<NUM_POINTS+1; i++)
    {
      scms = S_MIN + (S_MAX - S_MIN) * i / NUM_POINTS;
      mu_msbar = 2*pow(scms, 0.5);
      outFile << scms
              << "\t" << ru_pert_(scms, mu_msbar, alpha_s, mass_c, mass_b, mass_t)
              << "\t" << rd_pert_(scms, mu_msbar, alpha_s, mass_c, mass_b, mass_t)
              << "\t" << rs_pert_(scms, mu_msbar, alpha_s, mass_c, mass_b, mass_t)
              << "\t" << rc_pert_(scms, mu_msbar, alpha_s, mass_c, mass_b, mass_t) << std::endl;
    }

  outFile.close();

  return 0;
}
