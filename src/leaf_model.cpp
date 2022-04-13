#include <plant/leaf_model.h>

namespace plant {
Leaf::Leaf(): 
ci(1, NA_REAL),
g_c(1, NA_REAL),
A_lim(1, NA_REAL),
E(1, NA_REAL),
psi(1, NA_REAL),
profit(1, NA_REAL),
max_bound(2, NA_REAL),
min_bound(2, NA_REAL)
{}

double Leaf::calc_k_l_max(double K_s, double h_v, double h) const {
  return K_s * h_v / h;
}

double Leaf::calc_vul_b(double p_50, double c) const {
  return p_50/pow(-log(1 - 50.0/100.0),1/c);
}

// integrates
double Leaf::calc_cond_vuln(double psi, double k_l_max, double b,
                            double c) const {
  return k_l_max * exp(-pow((psi / b), c));
}

// replace f with some other function
double Leaf::calc_E_supply(double k_l_max, double b, double c, double psi_soil,
                           double psi_stem) {
  std::function<double(double)> f;
  f = [&](double psi) -> double { return calc_cond_vuln(psi, k_l_max, b, c); };

  return integrator.integrate(f, psi_soil, psi_stem);
}

double Leaf::calc_g_c(double psi_soil, double psi_stem, double k_l_max, double p_50,
                      double c, double b, double atm_kpa,
                      const double kg_2_mol_h2o, double atm_vpd) {
  double E_supply = calc_E_supply(k_l_max, c, b, psi_soil, psi_stem);

  return atm_kpa * E_supply * kg_2_mol_h2o / atm_vpd / 1.6;
}

double Leaf::calc_A_c(double ci_, double vcmax, double gamma_25,
                      double umol_per_mol_2_Pa, double km_25) {
  return (vcmax * (ci_ - gamma_25 * umol_per_mol_2_Pa)) / (ci_ + km_25);
}

double Leaf::calc_A_j(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                      double curv_fact, double a, double gamma_25,
                      double umol_per_mol_2_Pa, double ci_) {

  double jmax = vcmax * vcmax_25_2_jmax_25;
  double j = (a * PPFD + jmax -
              sqrt(pow(a * PPFD + jmax, 2) - 4 * curv_fact * a * PPFD * jmax)) /
             (2 * curv_fact);  // check brackets are correct

  return j / 4 *
         ((ci_ - gamma_25 * umol_per_mol_2_Pa) /
          (ci_ + 2 * gamma_25 * umol_per_mol_2_Pa));
}

double Leaf::calc_A_lim(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                        double curv_fact, double a, double gamma_25,
                        double umol_per_mol_2_Pa, double ci_, double km_25) {

  double A_c = calc_A_c(ci_, vcmax, gamma_25, umol_per_mol_2_Pa, km_25);
  double A_j = calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                        umol_per_mol_2_Pa, ci_);

  return (A_c + A_j - sqrt(pow(A_c + A_j, 2) - 4 * 0.98 * A_c * A_j)) / (2 * 0.98);
}

double Leaf::diff_ci(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                     double curv_fact, double a, double gamma_25,
                     double umol_per_mol_2_Pa, double x, double km_25,
                     double psi_soil, double psi_stem, double k_l_max, double p_50,
                     double c, double b, const double kg_2_mol_h2o,
                     double umol_per_mol_2_mol_per_mol, double atm_vpd,
                     double ca, double atm_kpa, double kPa_2_Pa) {

  double A_lim_ = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a,
                             gamma_25, umol_per_mol_2_Pa, x, km_25);

  double g_c_ = calc_g_c(psi_soil, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);

  return A_lim_ * umol_per_mol_2_mol_per_mol -
         (g_c_ * (ca - x) / (atm_kpa * kPa_2_Pa));
}

// need to fill in tol and max_iteratiosn
double Leaf::calc_assim_gross(double PPFD, double vcmax, double vcmax_25_2_jmax_25,
                         double curv_fact, double a, double gamma_25,
                         double umol_per_mol_2_Pa, double km_25,
                         double psi_soil, double psi_stem, double k_l_max,
                         double p_50, double c, double b,
                         const double kg_2_mol_h2o,
                         double umol_per_mol_2_mol_per_mol, double atm_vpd,
                         double ca, double atm_kpa, double kPa_2_Pa) {

  // not clear what x is here
  auto target = [&](double x) mutable -> double {
    return diff_ci(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                   umol_per_mol_2_Pa, x, km_25, psi_soil,
                   psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o,
                   umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa,
                   kPa_2_Pa);
  };

  // tol and iterations copied from control defaults (for now)
  ci[0] = util::uniroot(target, 0, 40, 1e-6, 1000);

  A_lim[0] = calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25,
                    umol_per_mol_2_Pa, ci[0], km_25);

  g_c[0] = calc_g_c(psi_soil, psi_stem, k_l_max, p_50, c,
                         b, atm_kpa, kg_2_mol_h2o, atm_vpd);

  E[0] = g_c[0] * 1.6 * atm_vpd / kg_2_mol_h2o / atm_kpa;

  psi[0] = psi_stem;                         

  return A_lim[0];
}


double Leaf::calc_hydraulic_cost(double psi_soil, double psi_stem, double k_l_max, double b, double c) {
  double k_l_soil_ = calc_cond_vuln(psi_soil, k_l_max, b, c);
  double k_l_stem_ = calc_cond_vuln(psi_stem, k_l_max, b, c);

  return k_l_soil_ - k_l_stem_;
}


double Leaf::calc_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
double psi_soil, double psi_stem, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, 
double umol_per_mol_2_mol_per_mol, double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
  double lambda_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_crit, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa)/
    calc_cond_vuln(psi_soil, k_l_max, b, c) - k_l_max * 0.05;
  double benefit_ = calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa);
  double cost_ = calc_hydraulic_cost(psi_soil, psi_stem, k_l_max, b, c);
  
  return benefit_ - lambda_*cost_;
}

// double Leaf::optimise_profit(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
//   double psi_soil, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, 
//   double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
//   psi[0] = psi_soil;

//   profit[0] = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi[0], k_l_max, 
//   p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 0) ;
    
//   if (psi_soil < psi_crit) {
//     double delta = 0.01;
//     double delta_crit = 1e-5;
//     psi[1] = psi[0] + delta;

//     profit[1] = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi[1], k_l_max, 
//       p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 1) ;

//     while (delta > delta_crit) {
      
//       psi[2] = psi[1] + delta;
//       profit[2]  = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi[2], k_l_max, 
//       p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 2) ;    
    
//       if ((profit[2] > profit[1]) && (profit[1] > profit[0])) {
//         profit[0] = profit[1];
//         profit[1] = profit[2];

//         psi[0] = psi[1];
//         psi[1] = psi[2];

//         A_lim[0] = A_lim[1];
//         A_lim[1] = A_lim[2];

//         g_c[0] = g_c[1];
//         g_c[1] = g_c[2];

//         ci[0] = ci[1];
//         ci[1] = ci[2];

//         E[0] = E[1];
//         E[1] = E[2];
      
//       } else {
//         delta = delta/2;
//       }
//     }
//   }

//   return(profit[0]);
//   }

double Leaf::optimise_profit_gss(double PPFD, double vcmax, double vcmax_25_2_jmax_25, double curv_fact, double a, double gamma_25, double umol_per_mol_2_Pa, double km_25, 
  double psi_soil, double k_l_max, double p_50, double c, double b, double kg_2_mol_h2o, double umol_per_mol_2_mol_per_mol, 
  double atm_vpd, double ca, double atm_kpa, double kPa_2_Pa, double psi_crit) {
  
  double gr = (sqrt(5) + 1)/2;

  //profit[0] = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, environment, psi_aep, b_CH, psi_soil_, k_l_max, 
  //p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, 0);

  min_bound[0] = psi_soil; 
  max_bound[0] = psi_crit;

  double opt_psi_stem = psi_soil;

  if (psi_soil < psi_crit) {
    double delta_crit = 1e-5;
    
    max_bound[1] = max_bound[0] - (max_bound[0] - min_bound[0])/gr;
    min_bound[1] = min_bound[0] + (max_bound[0] - min_bound[0])/gr;
  
    while (abs(max_bound[0] - min_bound[0]) > delta_crit) {
    
      if(calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, max_bound[1], k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit) >
      calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, min_bound[1], k_l_max, 
      p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit)) {
        max_bound[0] = min_bound[1];
      } else {
        min_bound[0] = max_bound[1];
      }

     min_bound[1] = min_bound[0] + (max_bound[0] - min_bound[0])/gr;
     max_bound[1] = max_bound[0] - (max_bound[0] - min_bound[0])/gr;


    }
    opt_psi_stem  = ((max_bound[0] + min_bound[0])/2);
    
  profit[0] = calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25,psi_soil, opt_psi_stem, k_l_max, 
  p_50, c, b, kg_2_mol_h2o, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit) ;

  }
  return(profit[0]);
}
}
// namespace plant
