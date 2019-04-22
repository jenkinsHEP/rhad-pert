#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "core/rhad_pert.h"
#include "core/Variable.h"

int main() {

/******************************************************************************/

/* Initialization */

  // Control parameters
  const int NUM_POINTS = 500;
  const int NUM_DICE = 1000;
  const double SCMS_MIN = 4.0;
  const double SCMS_MAX = 1000;
  /*
    NUM_POINTS     number of rows in output table
    NUM_DICE       number of monte carlo iterations
    SCMS_MIN       minimum center of mass energy squared in the scan
    SCMS_MAX       maximum center of mass energy squared in the scan
  */

  // Center of mass energy squared
  double scms;
  // msbar scale
  double mu_msbar;

  // Output stream
  std::ofstream outFile;

/******************************************************************************/

/* Parameters varied in error analysis */

  Variable1d *alpha_s_variable = new Variable1d(0.1181, 0.0011);
    double alpha_s_sample;
    double *alpha_s_histo = new double[NUM_DICE];

  Variable1d *mass_c_variable = new Variable1d(1.275, 0.025);
    double mass_c_sample;
    double *mass_c_histo = new double[NUM_DICE];

  Variable1d *mass_b_variable = new Variable1d(4.18, 0.04);
    double mass_b_sample;
    double *mass_b_histo = new double[NUM_DICE];

  Variable1d *mass_t_variable = new Variable1d(160.0, 5.0);
    double mass_t_sample;
    double *mass_t_histo = new double[NUM_DICE];

  Variable1d *log_mu_msbar_variable = new Variable1d(0.0, log(2.0));
    double mu_msbar_sample;
    double *log_mu_msbar_histo = new double[NUM_DICE];
  /*
    alpha_s		     strong coupling constant at the Z boson mass scale
    mass_c		     charm quark mass (msbar scale-independent)
    mass_b		     bottom quark mass (msbar scale-independent)
    mass_t		     top quark mass (msbar scale-independent)
    log_mu_msbar   variation of the logarithm of the msbar scale from its nominal value
  */

/* Results */
  double *ru_pert_histo = new double[NUM_DICE];
  double *rd_pert_histo = new double[NUM_DICE];
  double *rs_pert_histo = new double[NUM_DICE];
  double *rc_pert_histo = new double[NUM_DICE];

/******************************************************************************/

/* Print opening lines */

  outFile.open("rpert.txt");
  outFile << std::fixed << std::setprecision(6);

  outFile << "----------------------------------------------------------------------------------------------------------------------------------------\n\n"
          << " The contribution to the hadronic R-ratio from u,d,s,c quark flavors in perturbation theory. The columns should be read as follows:\n\n"
          << " s(GeV^2),  R_u,  d[R_u],  R_d,  d[R_d],  R_s,  d[R_s],  R_c,  d[R_c]\n\n"
          << " The following inputs were used (masses are in the MSbar scheme): \n\n"
          << " alpha_s(m_z)\t= " << alpha_s_variable->getMean() << " +/- " << alpha_s_variable->getError() << "\n"
          << " mass_c(mass_c)\t= " << mass_c_variable->getMean() << " +/- " << mass_c_variable->getError() << "\n"
          << " mass_b(mass_b)\t= " << mass_b_variable->getMean() << " +/- " << mass_b_variable->getError() << "\n"
          << " mass_t(mass_t)\t= " << mass_t_variable->getMean() << " +/- " << mass_t_variable->getError() << "\n\n"
          << " " << NUM_DICE << " iterations were used in the error analysis" << "\n"
          << " The MSbar scale was varied between 0.5 sqrt(s) and 2.0 sqrt(s)" << "\n\n"
          << "----------------------------------------------------------------------------------------------------------------------------------------\n\n";

/******************************************************************************/

/* Generate samples of parameters */

  for (int i=0; i<NUM_DICE; i++)
    {
      alpha_s_histo[i] = alpha_s_variable -> sample();
      mass_c_histo[i] = mass_c_variable -> sample();
      mass_b_histo[i] = mass_b_variable -> sample();
      mass_t_histo[i] = mass_t_variable -> sample();
      log_mu_msbar_histo[i] = log_mu_msbar_variable -> sample();
    }

  // Scan in center of mass energy squared
  for (int i=0; i<NUM_POINTS+1; i++)
    {
      // Set center of mass energy squared
      scms = SCMS_MIN + (SCMS_MAX - SCMS_MIN) * i / NUM_POINTS;

      // Calculate central values and uncertainties
      for (int j=0; j<NUM_DICE; j++)
        {
          mu_msbar = sqrt(scms) + exp( log_mu_msbar_histo[j] );

          ru_pert_histo[j] = ru_pert_(scms, mu_msbar, alpha_s_histo[j], mass_c_histo[j], mass_b_histo[j], mass_t_histo[j] );
          rd_pert_histo[j] = rd_pert_(scms, mu_msbar, alpha_s_histo[j], mass_c_histo[j], mass_b_histo[j], mass_t_histo[j] );
          rs_pert_histo[j] = rs_pert_(scms, mu_msbar, alpha_s_histo[j], mass_c_histo[j], mass_b_histo[j], mass_t_histo[j] );
          rc_pert_histo[j] = rc_pert_(scms, mu_msbar, alpha_s_histo[j], mass_c_histo[j], mass_b_histo[j], mass_t_histo[j] );
        }

      outFile << scms
              << "\t" << mean(ru_pert_histo, NUM_DICE) << "\t" << stdev(ru_pert_histo, NUM_DICE)
              << "\t" << mean(rd_pert_histo, NUM_DICE) << "\t" << stdev(rd_pert_histo, NUM_DICE)
              << "\t" << mean(rs_pert_histo, NUM_DICE) << "\t" << stdev(rs_pert_histo, NUM_DICE)
              << "\t" << mean(rc_pert_histo, NUM_DICE) << "\t" << stdev(rc_pert_histo, NUM_DICE) << std::endl;
    }

  outFile.close();

  return 0;
}
