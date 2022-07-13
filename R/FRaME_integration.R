#' Expands tidy_patch parameters
#' Copied from tidy_patch branch
#'
#' @param results Results from running tidy_patch
#'
#' @return
#' @export
#'

FF16_expand_state <- function(results) {
  
  data <- 
    results$species %>% 
    split(., .$species)
  
  for(i in seq_len(results$n_spp)) {
    
    s <- results$p$strategies[[i]]
    s$eta_c <- 1 - 2/(1 + s$eta) + 1/(1 + 2*s$eta)
    
    data[[i]] <- 
      data[[i]] %>%
      mutate(
        # These are formulas from ff16_strategy.cpp 
        # ideally wouldn't have to copy them here
        # could we expose them from startegy object
        # and call them directly?
        
        area_leaf = (height / s$a_l1)^(1.0 / s$a_l2),
        mass_leaf = area_leaf * s$lma,
        area_sapwood = area_leaf * s$theta,
        mass_sapwood = area_sapwood * height * s$eta_c * s$rho,
        area_bark = s$a_b1 * area_leaf * s$theta,
        mass_bark = area_bark * height * s$eta_c * s$rho,
        area_stem = area_bark + area_sapwood + area_heartwood,
        diameter_stem = sqrt(4 * area_stem / pi),
        mass_root = s$a_r1 * area_leaf,
        mass_live = mass_leaf + mass_sapwood + mass_bark + mass_root,
        mass_total =  mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root,
        mass_above_ground = mass_leaf + mass_bark + mass_sapwood +  mass_heartwood
      )
  }
  
  results$species <- data %>% bind_rows()
  
  results
}



#' Calculate crown parameters for a community at a specific age
#'
#' @param dat The results of run_scm_collect
#' @param tr A table of input traits
#' @param age Years since disturbance
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' 
#' @importFrom dplyr
#' @export
#'

shape_forest <- function(dat, tr, age, lat = -35, map = 1000, mat = 20) {
  
  result <- dat %>%
    plant::tidy_patch() %>% 
    plant::FF16_expand_state()
  
  # Run plantLitter model and add field
  ######################
  
  #Create function indiCrowns
  indiCrowns <- function(H, eta, lat = -35, map = 1000, mat = 20) {
    Hc <- H*0.91*(1-exp(-0.11*eta))
    He <- H*0.96*(1-exp(-0.21*eta))
    Ht <- H*0.96*(1-exp(-0.45*eta))
    W <- max(0.1*H, 0.2057*H+0.248*sqrt(abs(lat))+0.0005964*map+0.01202*mat-1.977)
    Af <- W*(Ht-He)+W*0.5*(H-Ht)+W*0.5*(He-Hc)
    
    res <- data.frame(matrix(ncol = 6, nrow = 1))
    x <- c("top", "base", "he", "ht", "W", "frontalArea")
    colnames(res) <- x
    res[1,]<-c(H, Hc, He, Ht, W, Af)
    
    return(res)
  }
  
  # Interpolate to set time
  full_data <- result$species
  study_age <- interpolate_to_times(full_data, age)
  study_age <- study_age[order(study_age$height),]%>%
    dplyr::filter(time == age, !is.na(height))
  
  # Add individual crown descriptors
  steps <- nrow(study_age)
  sp <- data.frame(matrix(ncol = 6, nrow = steps))
  colnames(sp) <- c("height", "base", "he", "ht", "w", "frontalArea")
  e <- as.vector(tr$eta) # Crown shape
  
  for (r in 1:steps) {
    sp[r,] <- indiCrowns(H=study_age$height[r], eta=e[as.numeric(study_age$species[r])],lat = lat, map = map, mat = mat)
  }
  all <- dplyr::left_join(study_age, sp, by = "height")
  
  # Update species names if tr table used
  #  if(!missing(tr)) {
  #    all <- left_join(all,tr, by = c("species" = "Species")) %>%
  #      mutate(species = name)
  #  }
  
  return(all)
}


#' Divide a plant community into strata
#' 
#' Creates a series of pseudo-transects to randomly sample the output data,
#' then divides these into 2-4 strata using k-means clustering, and
#' choosing the maximum number of significantly divided strata
#'
#' @param dat The results of run_scm_collect
#' @param tr A table of input traits
#' @param age Years since disturbance
#' @param lat Latitude (degrees)
#' @param map Mean annual precipitation (mm)
#' @param mat Mean annual temperature (degC)
#' @param propSamp Proportion of cohorts to test (0-1)
#' @param transects Number of repeats
#' 
#' @importFrom dplyr framer
#' @export
#'

stratify_community <- function(dat, tr, age, lat = -35, map = 1000, mat = 20, propSamp = 0.5, transects = 10){
  
  # Randomly sample results
  vegStrat <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(vegStrat) <- c("Species", "height", "base", "he", "ht")
  crowns <- plant:::shape_forest(dat, tr, age, lat, map, mat)
  steps <- nrow(crowns)
  samples <- ceiling(steps * propSamp)
  rand <- data.frame(matrix(ncol = transects, nrow = steps))
  for (t in 1:transects) {
    for (r in 1:steps) {
      # Randomise sampling
      rand[r,t] <- runif(n = 1, min = min(crowns$density), max = max(crowns$density))
    }
    # Order plants by the most likely to be sampled (cover x random)
    crowns$weight<-crowns$density*rand[,t]
    samp <- crowns[order(-crowns$weight),] %>%
      select(species, height, base, he, ht)
    # Select the chosen number of samples
    samp <- samp[1:samples, ]
    vegStrat <- rbind(vegStrat,samp)
  }
  
  # Calculate strata from random sample data
  vegStrat$Point <- 1
  vegStrat$top <- vegStrat$height
  vegStrat <- as.data.frame(vegStrat)
  vegStrat <- frameStratify(veg = vegStrat, spName = "species") %>%
    group_by(Stratum, species) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE)))%>%
    select(Stratum, species, height)
  #  vegStrat <- vegStrat[order(vegStrat$height),]
  
  # Interpolate to the stratum division heights, then add these to original data
  breaks <- nrow(vegStrat)-1
  heights <- vector()
  for (b in 1:breaks) {
    heights[b]<-vegStrat$height[b]
  }
  heights[breaks+1] <- max(crowns$height, na.rm = TRUE)
  heights <- as.numeric(heights)
  tidy_species_new <- interpolate_to_heights(crowns, heights) %>%
    select(!step)
  crowns <- crowns %>% select(!node)
  f<- rbind(tidy_species_new, crowns) %>% 
    drop_na() %>%
    mutate(Stratum = NA,
           d = NA,
           cohort = NA)
  
  # Measure density
  out <- f[0,]
  outA <- f[0,]
  for (sp in unique(vegStrat$species)) {
    s <- filter(vegStrat, species == sp)
    outA <- filter(f, species == sp)
    s$height[nrow(s)] <- 1000
    # Identify the stratum of each plant
    outA$Stratum <- cut(as.numeric(outA$height), breaks = c(0,s$height),
                        labels = s$Stratum,
                        include.lowest = TRUE)
    # Integrate density
    outA <- outA[order(outA$height),]
    for (h in 1:nrow(outA)) {
      if (h == 1) {
        outA$d[h] <- max(((outA$density[h]/2) * outA$height[h]),0)
      } else {
        outA$d[h] <- max((mean(outA$density[h], outA$density[h-1]) * (outA$height[h] - outA$height[h-1])),0)
      }
    }
    
    out <- rbind(out, outA)
    out <- out[which(out$d>0),] #Remove any empty cohorts
  }
  return(out)
}

#' Wrapper to grow a forest in plant ready for fire modelling
#' 
#' Grows the forest from an imported table of traits
#' Runs tidy_patch and expand_state
#'
#' @param dat An input table listing species traits
#' @param B_lf1 Potential CO2 photosynthesis at average leaf nitrogen
#' @param dist Mean disturbance interval for investigating trends
#' 
#' @importFrom dplyr
#' @export
#'
grow_forest <- function(dat, B_lf1 = 0.8273474, dist = 1000) {
  
  plant::plant_log_console()
  
  # Collect plant traits
  params <- scm_base_parameters("FF16")
  #  param$birth_rate <- 
  params$disturbance_mean_interval <- dist
  eta <- as.vector(dat$eta) # Crown shape
  lma <- as.vector(dat$lma) # kgm−2
  rho <- as.vector(dat$rho) # Wood density - kgm−3
  theta <- as.vector(dat$theta)  # Sapwood per unit leaf area
  a_l1 <- as.vector(dat$a_l1) # Height of plant with leaf area of 1m2
  a_l2 <- as.vector(dat$a_l2) # Exponent of relationship between height and leaf area
  a_b1 <- as.vector(dat$a_b1) # Ratio of bark area to sapwood area
  hmat <- as.vector(dat$hmat) # Height at maturation (m)
  a_f1 <- as.vector(dat$a_f1) # Maximum allocation to reproduction
  a_f2 <- as.vector(dat$a_f2) # Parameter determining rate of change in r(x,ml) around Hmat
  FF16_hyperpar1 <- make_FF16_hyperpar(B_lf1 = B_lf1) # Potential CO2 photosynthesis at average leaf Nitrogen - mold−1 m−2
  
  patch <- expand_parameters(trait_matrix(x = c(eta, lma, rho, theta, a_b1, hmat, a_l1, a_l2, a_f1, a_f2), 
                                          trait_name = c("eta", "lma", "rho", "theta", "a_b1", "hmat", "a_l1", "a_l2", "a_f1", "a_f2")), 
                             params, 
                             FF16_hyperpar1,
                             birth_rate_list = dat$birth_rate,
                             mutant = FALSE)
  
  patch <- build_schedule(patch)
  dat <- run_scm_collect(patch)
  
  return(dat)
}


#' Builds randomised forests for fire modelling
#' @export

randomise_forest <- function(){
  nSp <- round(runif(1, min = 1, max = 4), 0) # Number of species in the forest, between 1 & 4
  
  # Build the table
  spTab <- data.frame(matrix(ncol = 13, nrow = nSp))
  colnames(spTab) <- c("name", "Species", "eta", "lma", "rho", "a_b1", "a_l1", "a_l2", "a_f1", "a_f2", "hmat", "theta", "birth_rate")
  
  spTab$name <- letters[1:nSp]                          # Name for the species
  spTab$Species <- as.character(seq(1:nSp))             # Number of the species
  spTab$eta <- runif(nSp, min = 1, max = 12)            # Crown shape parameter
  spTab$lma <- runif(nSp, min = 0.007, max = 2.6)       # Leaf mass per area
  spTab$rho <- runif(nSp, min = 200, max = 1100)        # Wood density
  spTab$a_b1 <- runif(nSp, min = 0.085, max = 0.34)     # Ratio of bark area to sapwood area
  spTab$a_l1 <- runif(nSp, min = 0, max = 40)           # Height of plant with leaf area of 1m2
  spTab$a_l2 <- runif(nSp, min = 0.2, max = 0.4)        # Exponent of relationship between height and leaf area
  spTab$a_f1 <- runif(nSp, min = 0, max = 1)            # Maximum allocation to reproduction
  spTab$a_f2 <- runif(nSp, min = 0, max = 50)           # Rate of change in reproduction
  spTab$hmat <- runif(nSp, min = 0.1, max = 50)         # Height at maturation
  spTab$theta <- runif(nSp, min = 0.0001, max = 0.00025) # Sapwood area per unit leaf area
  spTab$leaf_thickness <- NA
  spTab$leafForm <- NA
  spTab$moisture <- NA
  spTab$propDead <- NA
  spTab$ignitionTemp <- NA
  spTab$lwRat <- NA
  spTab$leaf_area <- NA
  spTab$bark_density <- NA
  spTab$G.C_rat <- NA
  spTab$C.C_rat <- NA
  spTab$stemOrder <- NA
  spTab$birth_rate <- 1
  
  return(spTab)
}