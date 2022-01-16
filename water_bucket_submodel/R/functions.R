calc_soil_wc <- function(theta, theta_sat){
  soil_wc <- theta/theta_sat
  soil_wc
}

calc_psi <- function(a, n, soil_wc) {
  psi <- -a*(soil_wc)^(-n)
  psi
}

calc_emmissivity_air_Kang <- function(air_temp){
  emmissivity_air_Kang <- 1.24*(2/air_temp)^(1/7) # Comparison and analysis of bare soil evaporation models combined with ASTER data in Heihe River Basin
  emmissivity_air_Kang
}

calc_diffus_water_vapour <- function(soil_temp) {
  diffus_water_vapour <- 24.2*(10^-6)*(soil_temp/293.2)^1.75 #m^2 s^-1 Jones 1992 appendix
  diffus_water_vapour
  }

calc_eff_diffus_water_vapour <- function(diffus_water_vapour, theta_sat, tortuosity) {
  eff_diffus_water_vapour <- diffus_water_vapour*theta_sat*tortuosity
  eff_diffus_water_vapour
  }

calc_dry_layer <- function(soil_wc, depth_top_layer) {
  dry_layer<-depth_of_top_layer*(1-soil_wc)
  dry_layer
}

calc_G_ws <- function(eff_diffus_water_vapour, theta_sat, dry_layer) {
  G_ws <- eff_diffus_water_vapour*(theta_sat/dry_layer)
  G_ws
}

calc_G_am <- function(u,z,d,zo){
G_am <- von_Karman^2*u/(log((z-d)/zo))^2 # Jones 1992, will need to get wind profile for a forest which will supply values, also check MAESPA, this is the atmospheric conductance for the ground
G_am
}

calc_E_sat <- function(temp){
  temp_C <- temp-273.15
  a <- 611.21
  b <- 18.678 - (temp_C/234.5)
  c <- 257.14
  f <- (1.0007+3.46*10^-8)
  e_sat_pa <- f*(a*exp(b*temp_C/(c+temp_C)))
  e_sat_pa
}

calc_E_soil <- function(psi, soil_temp, E_sat){
  E_sat*kPa_to_Pa*exp(psi*Pa_to_MPa*partial_molal_water_vapour/(molar_gas_constant*soil_temp))
}

calc_rho <- function(temp){
  rho <- 353/temp
  rho
}

