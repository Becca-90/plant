// -*-c++-*-
#ifndef PLANT_PLANT_LEAF_MODEL_H_
#define PLANT_PLANT_LEAF_MODEL_H_

// TODO: replace with constants
// #define umol_per_mol_2_mol_per_mol 0.000...
// #define umol_per_mol_2_Pa ...
// #define kg_2_mol_h2o ...
// #define kPa_2_Pa ...

#include <plant/models/ff16_environment.h>
#include <plant/qag.h>
#include <plant/uniroot.h>

namespace plant {

class Leaf {
public:
  Leaf();

  quadrature::QAG integrator;

  double ci;
  double g_c;
  double A_lim;
  double E;
  double psi;
  double profit;





  // this might end up hard-coded
  void initialize_integrator(int integration_rule = 21, double integration_tol = 1e-6) {

    integrator = quadrature::QAG(integration_rule,
                                 1, // fixed integration
                                 integration_tol, integration_tol);
  }

  double calc_k_l_max(double K_s, double h_v, double h) const;

  double calc_vul_b(double p_50, double c) const;
  double calc_cond_vuln(double psi, double k_l_max, double b, double c) const;

  double calc_E_supply(double k_l_max, double b, double c, double psi_soil,
                       double psi_stem);

  double calc_g_c(double psi_soil, double psi_stem, double k_l_max, double p_50,
                  double c, double b, double atm_kpa, const double kg_2_mol_h2o,
                  double atm_vpd); // define as a constant

  double calc_A_c(double ci_, double vcmax, double gamma_25,
                  double umol_per_mol_2_Pa, double km_25);
  double calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                  double curv_fact, double a, double gamma_25,
                  double umol_per_mol_2_Pa, double ci_);
  double calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                    double curv_fact, double a, double gamma_25,
                    double umol_per_mol_2_Pa, double ci_, double km_25);

  double diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                     double curv_fact, double a, double gamma_25,
                     double umol_per_mol_2_Pa, double x, double km_25,
                     double psi_soil, double psi_stem, double k_l_max, double p_50,
                     double c, double b, const double kg_2_mol_h2o,
                     double umol_per_mol_2_mol_per_mol, double atm_vpd,
                     double ca, double atm_kpa, double kPa_2_Pa);
    
  double calc_assim_gross(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         double psi_soil, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa);
  double calc_hydraulic_cost(double psi_soil, double psi_stem, double k_l_max, double b, double c);
  double calc_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25, double psi_soil, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit);                                                 
double optimise_profit_gss(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         double psi_soil, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit);
double calc_hydraulic_cost_bartlett(double psi_soil, double psi_stem, double k_l_max, double b, double c, double beta, double beta_2, double huber_value, double height);
  double calc_profit_bartlett(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25, double psi_soil, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit, double beta, double beta_2, double huber_value, double height);                                                 
double optimise_profit_gss_bartlett(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         double psi_soil, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa, double psi_crit, double beta, double beta_2, double huber_value, double height);                          
                          
};

} // namespace plant
#endif