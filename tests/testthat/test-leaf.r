source("tests/reference_leaf.R")
context("SCM-general")


test_that("Basic functions", {
  
  l <- Leaf()
  
  library(tidyverse)
  
  K_s = 2 #stem-specific conductivity (kg h2o m^-1 stem s^-1 MPa^-1)
  h_v = 0.000157 #huber value (m^2 sapwood area m^-2 leaf area)
  h = 5 #height/path length (m)
  
  expect_equal(l$calc_k_l_max(K_s, h_v, h), calc_k_l_max(K_s, h_v, h))
  
  p_50 = 1.731347 #stem water potential at 50% loss of conductivity
  c = 2.04 #shape parameter for hydraulic vulnerability curve (unitless) estimated from trait data in Austraits from Choat et al. 2012
  
  expect_equal(l$calc_vul_b(p_50, c), calc_vul_b(p_50, c))
  
  k_l_max <- l$calc_k_l_max(K_s, h_v, h)  #leaf specific conductance value (kg m^-2 leaf area m^-1 s^-1 MPa^-1)
  psi <- 1 #nominated value for water potential for testing vulnerability curve equations only (-MPa)
  b <- 2.072101 #shape parameter for vulnerability curve, point of 37% conductance (-MPa) 
  
  expect_equal(l$calc_cond_vuln(psi, k_l_max, b, c), calc_k_l(psi, k_l_max, b, c), tolerance = 1e-8)
  
  psi_soil <- 0.5 #nominated value for soil water potential for testing (-MPa)
  
  atm_vpd <- 1 #atmopsheric vapour pressure deficit (kPa)
  psi_stem <- psi_soil+1 #stem water potential [set at above soil water potential for testing] (-MPa)
  
  #for situations where psi_soil is < than psi_crit

  expect_equal(l$calc_E_supply(k_l_max, b, c, psi_soil, psi_stem), calc_E_supply(psi_stem, psi_soil, k_l_max = k_l_max, b=b, c=c)$value)
  expect_equal(l$calc_g_c(psi_soil, psi_stem, k_l_max, p_50, c, b, atm_kpa, kg_2_mol_h20, atm_vpd), calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b))
  
  vcmax = 100 #maximum carboxylation rate, defined by leaf nitrogen (umol m^-2 s^-1) 
  c_i = 30 #intra-cellular carbon dioxide parital pressure (Pa)

  expect_equal(l$calc_A_c(c_i, vcmax, gamma_25, umol_per_mol_2_Pa, km_25), (vcmax * (c_i - gamma_25*umol_per_mol_2_Pa))/(c_i + km_25))

  PPFD <- 900 #Photon flux density (umol m^-2 s^-1)
  ca <- 40 #atmospheric carbon dioxide partial pressure
  
  psi_crit <- calc_psi_crit(b, c) #stem water potential at which conductance is 95%
  
  expect_equal(l$calc_A_c(c_i, vcmax, gamma_25, umol_per_mol_2_Pa, km_25), calc_A_c(c_i = c_i, vcmax = vcmax))
  
  expect_equal(l$calc_A_j(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i), calc_A_j(c_i = c_i, PPFD = PPFD, vcmax = vcmax))
  
  expect_equal(l$calc_A_lim(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, c_i, km_25), calc_A_lim(c_i, vcmax, PPFD))
  
  expect_equal(l$calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa), calc_ben_gross(psi_soil = psi_soil, psi_stem = psi_stem,  b = b, c = c, k_l_max = k_l_max, PPFD = PPFD, vcmax = vcmax, ca =ca, atm_vpd = atm_vpd), tolerance = 1e-6)
  
  expect_equal(l$calc_hydraulic_cost(psi_soil, psi_stem, k_l_max, b, c), calc_hydraulic_cost(psi_soil = psi_soil, psi_stem = psi_stem,  b = b, c = c, k_l_max = k_l_max))
  
  expect_equal(l$calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h20, 
                umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),calc_profit(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca))
  
  expect_equal(l$E,calc_E_supply(psi_stem, psi_soil, k_l_max = k_l_max, b=b, c=c)$value)
  
  expect_equal(l$ci, solve_ci(vcmax, PPFD, calc_g_c(psi_stem = psi_stem, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b), ca)[[1]])
  
  expect_equal(l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                     psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                     atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),
               find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$profit_root)
  
  
  expect_equal(l$E,find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$E_root)
  
  expect_equal(l$ci,find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$ci_root)
  
  beta = 15000 #cost term for Bartlett profit equation (umol m^-3 stem s^-1)
  beta_2 = 1 #shape term for Bartlett profit equation
  
  expect_equal(l$calc_hydraulic_cost_bartlett(psi_soil, psi_stem, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h),calc_hydraulic_cost_bartlett(psi_soil, psi_stem, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h))
  
  expect_equal(l$calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                      psi_soil, psi_stem, k_l_max, p_50, c, b, kg_2_mol_h20, 
                                      umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value = h_v, height = h),
               calc_profit_bartlett(psi_stem, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h))

  expect_equal(l$optimise_profit_gss_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                 psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                 atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value = h_v, height = h),
  find_opt_psi_stem_bartlett(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h)$profit_root)
  
  ##### when psi soil is equal to psi stem
  
  
  expect_equal(l$calc_E_supply(k_l_max, b, c, psi_soil, psi_soil), calc_E_supply(psi_soil, psi_soil, k_l_max = k_l_max, b=b, c=c)$value)
  
  expect_equal(l$calc_g_c(psi_soil, psi_soil, k_l_max, p_50, c, b, atm_kpa, kg_2_mol_h20, atm_vpd), calc_g_c(psi_stem = psi_soil, psi_soil = psi_soil, atm_vpd = atm_vpd, k_l_max = k_l_max, c = c, b = b))
  
  expect_equal(l$calc_assim_gross(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                  psi_soil, psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, 
                                  umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa),calc_ben_gross(psi_soil = psi_soil, psi_stem = psi_soil, atm_vpd= atm_vpd, k_l_max =k_l_max, b = b, c = c, PPFD = PPFD, vcmax = vcmax, ca = ca))
  
  expect_equal(calc_lambda(psi_crit, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca),0)
  
  expect_equal(l$calc_profit(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                             psi_soil, psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, 
                             umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),calc_profit(psi_soil, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca))
  
  expect_equal(l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                     psi_crit, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                     atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),find_opt_psi_stem(psi_soil = psi_crit, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$profit_root)
  
  expect_equal(l$ci,find_opt_psi_stem(psi_soil = psi_crit, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$ci_root)
  
  
  
  expect_equal(l$calc_hydraulic_cost_bartlett(psi_soil, psi_soil, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h),calc_hydraulic_cost_bartlett(psi_soil, psi_soil, k_l_max, b, c, beta, beta_2, huber_value = h_v, height = h))
  
  expect_equal(l$calc_profit_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                      psi_soil, psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, 
                                      umol_per_mol_2_mol_per_mol, atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value = h_v, height = h),
               calc_profit_bartlett(psi_soil, psi_soil, k_l_max, vcmax, PPFD, b, c, psi_crit, atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h))
  
  expect_equal(l$optimise_profit_gss_bartlett(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                              psi_crit, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                              atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit, beta, beta_2, huber_value = h_v, height = h),
               find_opt_psi_stem_bartlett(psi_soil = psi_crit, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca, beta, beta_2, huber_value = h_v, height = h)$profit_root)
  
  ###### when psi soil is greater than psi crit
  
  psi_soil = psi_crit + 1
  
  expect_equal(l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                    psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                    atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD, b = b, c = c, psi_crit = psi_crit, atm_vpd = atm_vpd, ca = ca)$profit_root)
  
  expect_equal(l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                     psi_crit, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                     atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit),
               l$optimise_profit_gss(PPFD, vcmax, vcmax_25_2_jmax_25, curv_fact, a, gamma_25, umol_per_mol_2_Pa, km_25, 
                                     psi_soil, k_l_max, p_50, c, b, kg_2_mol_h20, umol_per_mol_2_mol_per_mol, 
                                     atm_vpd, ca, atm_kpa, kPa_2_Pa, psi_crit))
  
  }
)
