# CUSTOM FUNCTIONS

# Effect size (lnVR and lnCVR) for 2 main effects and interaction effect ----
effect_setV <- function(CC_n, CC_mean, CC_SD,
                       EC_n, EC_mean, EC_SD,
                       CS_n, CS_mean, CS_SD,
                       ES_n, ES_mean, ES_SD,
                       percent){
  
  if(percent == "no"){
  # lnRR----
  # main effect Environmental enrichment----
  lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                         log(0.5*(CS_mean+ CC_mean))
  
  lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
    (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
  
  # main effect Stress----
  lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                         log(0.5*(EC_mean+ CC_mean))
  
  lnRRV_S <- lnRRV_E
  
  # interaction----
  
  lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                            (log(EC_mean) - log(CC_mean))
  
  
  lnRRV_ES <- 
    (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
     ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
      ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
       ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
  
  # SMD
  SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                                (EC_n-1)*EC_SD^2 + 
                                (CS_n-1)*CS_SD^2 +
                                (CC_n-1)*CC_SD^2) / 
                               (ES_n + EC_n + CS_n + CC_n - 4))
  
  
  
  # lnVR
  # main effect Environmental enrichment----
  lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_S <- lnVRV_E
  
  # interaction----
  
  lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV <- ES_SD/ES_mean
  EC_CV <- EC_SD/EC_mean
  CS_CV <- CS_SD/CS_mean
  CC_CV <- CS_SD/CS_mean
  
  lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_E <- lnRRV_E + lnVRV_E 
  
  # main effect Stress----
  lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_S <- lnRRV_S + lnVRV_S
  
  # interaction----
  lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
  effect <- tibble(
    # lnVR
    lnVR_E = lnVR_E,
    lnVRV_E = lnVRV_E, 
    lnVR_S = lnVR_S, 
    lnVRV_S = lnVRV_S,
    lnVR_ES =lnVR_ES, 
    lnVRV_ES = lnVRV_ES,
    # lnCVR
    lnCVR_E = lnCVR_E,
    lnCVRV_E = lnCVRV_E, 
    lnCVR_S = lnCVR_S, 
    lnCVRV_S = lnCVRV_S,
    lnCVR_ES =lnCVR_ES, 
    lnCVRV_ES = lnCVRV_ES
  )
  effect
  }
  
  else {
    
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
    
    # lnRR
    # main effect Enrichment
    lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                           log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) +  
                             (1/(CS_mean + CC_mean))^2*(CS_SD^2 /CS_n + CC_SD^2 / CC_n) 
    
    # main effect Stress
    lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                           log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_S <- lnRRV_E
    
    # interaction----
    
    lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                              (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ES <- (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
                    ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
                    ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
                    ((CC_SD)^2 / ((CC_mean)^2*CC_n)))    
     
     
    # lnVR
    # main effect Environmental enrichment----              
    lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
    # main effect Stress----
    lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnVRV_S <- lnVRV_E
  
    # interaction----
  
    lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    # lnCVR
    # main effect Environmental enrichment----
    ES_CV = ES_SD/ES_mean
    EC_CV = EC_SD/EC_mean
    CS_CV = CS_SD/CS_mean
    CC_CV = CS_SD/CS_mean
  
    lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnCVRV_E <- lnRRV_E + lnVRV_E 
  
    # main effect Stress----
    lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnCVRV_S <- lnRRV_S + lnVRV_S
  
    # interaction----
  
    lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
    effect <- tibble(
    # lnVR
    lnVR_E = lnVR_E,
    lnVRV_E = lnVRV_E, 
    lnVR_S = lnVR_S, 
    lnVRV_S = lnVRV_S,
    lnVR_ES =lnVR_ES, 
    lnVRV_ES = lnVRV_ES,
    # lnCVR
    lnCVR_E = lnCVR_E,
    lnCVRV_E = lnCVRV_E, 
    lnCVR_S = lnCVR_S, 
    lnCVRV_S = lnCVRV_S,
    lnCVR_ES = lnCVR_ES, 
    lnCVRV_ES = lnCVRV_ES
  )
    effect
  }
  
}



# Effect size (lnRR and SMD) for 2 main effects and interaction effect ----
effect_set <- function(CC_n, CC_mean, CC_SD,
                       EC_n, EC_mean, EC_SD,
                       CS_n, CS_mean, CS_SD,
                       ES_n, ES_mean, ES_SD,
                       percent){
  
  if(percent == "no"){
  
  # lnRR----
  # main effect Environmental enrichment----
  lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                         log(0.5*(CS_mean+ CC_mean))
  
  lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
    (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
  
  # main effect Stress----
  lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                         log(0.5*(EC_mean+ CC_mean))
  
  lnRRV_S <- lnRRV_E
  
  # interaction----
  
  lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                            (log(EC_mean) - log(CC_mean))
  
  
  lnRRV_ES <- 
    (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
     ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
      ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
       ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
  
  # SMD
  SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                                (EC_n-1)*EC_SD^2 + 
                                (CS_n-1)*CS_SD^2 +
                                (CC_n-1)*CC_SD^2) / 
                               (ES_n + EC_n + CS_n + CC_n - 4))
  
  
  # main effect Environment enrichment
  SMD_E <- ((ES_mean + EC_mean) - (CS_mean + CC_mean))/ (2*SD_pool)
  
  
  SMDV_E <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + 
                    (SMD_E^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
  
  
  
  # main effect Stress
  SMD_S <- ((ES_mean + CS_mean) - (EC_mean + CC_mean)) / (2*SD_pool)
  
  SMDV_S <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + 
                    (SMD_S^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
  
  # interaction
  SMD_ES <- ((ES_mean - EC_mean) - (CS_mean - CC_mean)) / SD_pool
  
  SMDV_ES <- (1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_ES^2 / (2*(ES_n + EC_n + CS_n + CC_n)))
  
  
  # lnVR
  # main effect Environmental enrichment----
  lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_S <- lnVRV_E
  
  # interaction----
  
  lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV <- ES_SD/ES_mean
  EC_CV <- EC_SD/EC_mean
  CS_CV <- CS_SD/CS_mean
  CC_CV <- CS_SD/CS_mean
  
  lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_E <- lnRRV_E + lnVRV_E 
  
  # main effect Stress----
  lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_S <- lnRRV_S + lnVRV_S
  
  # interaction----
  lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
  effect <- tibble(
    # lnRR
    lnRR_E = lnRR_E,
    lnRRV_E = lnRRV_E, 
    lnRR_S = lnRR_S, 
    lnRRV_S = lnRRV_S,
    lnRR_ES =lnRR_ES, 
    lnRRV_ES = lnRRV_ES,
    #SMD
    SMD_E = SMD_E,
    SMDV_E = SMDV_E, 
    SMD_S = SMD_S, 
    SMDV_S = SMDV_S, 
    SMD_ES = SMD_ES, 
    SMDV_ES = SMDV_ES
    # lnVR
    #lnVR_E = lnVR_E,
    #lnVRV_E = lnVRV_E, 
    #lnVR_S = lnVR_S, 
    #lnVRV_S = lnVRV_S,
    #lnVR_ES =lnVR_ES, 
    #lnVRV_ES = lnVRV_ES,
    # lnCVR
    #lnCVR_E = lnCVR_E,
    #lnCVRV_E = lnCVRV_E, 
    #lnCVR_S = lnCVR_S, 
    #lnCVRV_S = lnCVRV_S,
    #lnCVR_ES =lnCVR_ES, 
    #lnCVRV_ES = lnCVRV_ES
  )
  effect
  }
  
  else {
    
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
     
    # lnRR
    # main effect Enrichment
    lnRR_E <- log(0.5*(ES_mean + EC_mean)) - 
                           log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_E <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) +  
                             (1/(CS_mean + CC_mean))^2*(CS_SD^2 /CS_n + CC_SD^2 / CC_n) 
    
    # main effect Stress
    lnRR_S <- log(0.5*(ES_mean + CS_mean)) - 
                           log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_S <- lnRRV_E
    
    # interaction----
    
    lnRR_ES <-   (log(ES_mean) - log(CS_mean)) - 
                              (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ES <- (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
                    ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
                    ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
                    ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
    
    # SMD
    SD_pool <- sqrt(((ES_n-1)*ES_SD^2 + 
                    (EC_n-1)*EC_SD^2 + 
                    (CS_n-1)*CS_SD^2 +
                    (CC_n-1)*CC_SD^2) / 
                    (ES_n + EC_n + CS_n + CC_n - 4))
    
    
    # main effect Environment enrichment
    SMD_E <- ((ES_mean + EC_mean) - (CS_mean + CC_mean))/ (2*SD_pool)
    
    
    SMDV_E <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_E^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
    
    # main effect Stress
    SMD_S <- ((ES_mean + CS_mean) - (EC_mean + CC_mean)) / (2*SD_pool)
    
    SMDV_S <- 0.25*((1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_S^2 /(2*(ES_n + EC_n + CS_n + CC_n))))
    
    # interaction
    SMD_ES <- ((ES_mean - EC_mean) - (CS_mean - CC_mean)) / SD_pool
    
    SMDV_ES <- (1/ES_n + 1/EC_n + 1/CS_n + 1/CC_n) + (SMD_ES^2 / (2*(ES_n + EC_n + CS_n + CC_n)))
    
    # lnVR
    # main effect Environmental enrichment----              
    lnVR_E <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnVRV_E <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
    # main effect Stress----
    lnVR_S <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnVRV_S <- lnVRV_E
  
    # interaction----
  
    lnVR_ES <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnVRV_ES <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    # lnCVR
    # main effect Environmental enrichment----
    ES_CV = ES_SD/ES_mean
    EC_CV = EC_SD/EC_mean
    CS_CV = CS_SD/CS_mean
    CC_CV = CS_SD/CS_mean
  
    lnCVR_E <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
    lnCVRV_E <- lnRRV_E + lnVRV_E 
  
    # main effect Stress----
    lnCVR_S <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
    lnCVRV_S <- lnRRV_S + lnVRV_S
  
    # interaction----
  
    lnCVR_ES <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
    lnCVRV_ES <- lnRRV_ES + lnVRV_ES
  
    effect <- tibble(
    # lnRR
    lnRR_E = lnRR_E,
    lnRRV_E = lnRRV_E, 
    lnRR_S = lnRR_S, 
    lnRRV_S = lnRRV_S,
    lnRR_ES =lnRR_ES, 
    lnRRV_ES = lnRRV_ES,
    #SMD
    SMD_E = SMD_E,
    SMDV_E = SMDV_E, 
    SMD_S = SMD_S, 
    SMDV_S = SMDV_S, 
    SMD_ES = SMD_ES, 
    SMDV_ES = SMDV_ES
    # lnVR
    #lnVR_E = lnVR_E,
    #lnVRV_E = lnVRV_E, 
    #lnVR_S = lnVR_S, 
    #lnVRV_S = lnVRV_S,
    #lnVR_ES =lnVR_ES, 
    #lnVRV_ES = lnVRV_ES,
    # lnCVR
    #lnCVR_E = lnCVR_E,
    #lnCVRV_E = lnCVRV_E, 
    #lnCVR_S = lnCVR_S, 
    #lnCVRV_S = lnCVRV_S,
    #lnCVR_ES = lnCVR_ES, 
    #lnCVRV_ES = lnCVRV_ES
  )
    effect
  }
  
}


# Removing asin_trans for sensitivity analysis----

effect_setb <- function(CC_n, CC_mean, CC_SD,
                        EC_n, EC_mean, EC_SD,
                        CS_n, CS_mean, CS_SD,
                        ES_n, ES_mean, ES_SD)
  {
    # lnRR----
    # main effect Environmental enrichment----
    lnRR_Eb <- log(0.5*(ES_mean + EC_mean)) - 
      log(0.5*(CS_mean+ CC_mean))
    
    lnRRV_Eb <-  (1/(ES_mean + EC_mean))^2*(ES_SD^2 / ES_n + EC_SD^2 / EC_n) + 
      (1/(CS_mean + CC_mean))^2*(CS_SD^2 / CS_n + CC_SD^2 / CC_n)
    
    # main effect Stress----
    lnRR_Sb <- log(0.5*(ES_mean + CS_mean)) - 
      log(0.5*(EC_mean+ CC_mean))
    
    lnRRV_Sb <- lnRRV_Eb
    
    # interaction----
    
    lnRR_ESb <-   (log(ES_mean) - log(CS_mean)) - 
      (log(EC_mean) - log(CC_mean))
    
    
    lnRRV_ESb <- 
      (((ES_SD)^2 / ((ES_mean)^2*ES_n)) + 
         ((EC_SD)^2 / ((EC_mean)^2*EC_n)) + 
         ((CS_SD)^2 / ((CS_mean)^2*CS_n)) +
         ((CC_SD)^2 / ((CC_mean)^2*CC_n)))
         
         
  # lnVR
  # main effect Environmental enrichment----
                         
  lnVR_Eb <- 0.5*log((ES_SD*EC_SD)/(CS_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnVRV_Eb <- 0.25*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1))
  
  # main effect Stress----
  lnVR_Sb <- 0.5*log((ES_SD*CS_SD)/(EC_SD*CC_SD)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnVRV_Sb <- lnVRV_Eb
  
  # interaction----
  
  lnVR_ESb <- log((ES_SD/CS_SD)/(EC_SD/CC_SD)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnVRV_ESb <- 0.5/(ES_n - 1) + 0.5/(EC_n - 1) + 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  # lnCVR
  # main effect Environmental enrichment----
  ES_CV = ES_SD/ES_mean
  EC_CV = EC_SD/EC_mean
  CS_CV = CS_SD/CS_mean
  CC_CV = CS_SD/CS_mean
  
  lnCVR_Eb <- 0.5*log((ES_CV*EC_CV)/(CS_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) + 0.5/(EC_n - 1) - 0.5/(CS_n - 1) - 0.5/(CC_n - 1))                     
  
  lnCVRV_Eb <- lnRRV_Eb + lnVRV_Eb 
  
  # main effect Stress----
  lnCVR_Sb <- 0.5*log((ES_CV*CS_CV)/(EC_CV*CC_CV)) + 0.5*(0.5/(ES_n - 1) - 0.5/(EC_n - 1) + 0.5/(CS_n - 1) - 0.5/(CC_n - 1))
  
  lnCVRV_Sb <- lnRRV_Sb + lnVRV_Sb
  
  # interaction----
  
  lnCVR_ESb <- log((ES_CV/CS_CV)/(EC_CV/CC_CV)) + 0.5/(ES_n - 1) - 0.5/(EC_n - 1) - 0.5/(CS_n - 1) + 0.5/(CC_n - 1)
  
  lnCVRV_ESb <- lnRRV_ESb + lnVRV_ESb
  
        
    effectb <- tibble(
      # lnRR
      lnRR_Eb = lnRR_Eb,
      lnRRV_Eb = lnRRV_Eb, 
      lnRR_Sb = lnRR_Sb, 
      lnRRV_Sb = lnRRV_Sb,
      lnRR_ESb =lnRR_ESb, 
      lnRRV_ESb = lnRRV_ESb,
      # lnVR
      lnVR_Eb = lnVR_Eb,
      lnVRV_Eb = lnVRV_Eb, 
      lnVR_Sb = lnVR_Sb, 
      lnVRV_Sb = lnVRV_Sb,
      lnVR_ESb =lnVR_ESb, 
      lnVRV_ESb = lnVRV_ESb,
      # lnCVR
      lnCVR_Eb = lnCVR_Eb,
      lnCVRV_Eb = lnCVRV_Eb, 
      lnCVR_Sb = lnCVR_Sb, 
      lnCVRV_Sb = lnCVRV_Sb,
      lnCVR_ESb =lnCVR_ESb, 
      lnCVRV_ESb = lnCVRV_ESb
    )
    effectb
}


# Pairwise comparisons lnRR (not for SMD) -----

effect_set2 <- function(CC_n, CC_mean, CC_SD,
                        EC_n, EC_mean, EC_SD,
                        CS_n, CS_mean, CS_SD,
                        ES_n, ES_mean, ES_SD,
                        percent){
  
  if(percent == "no"){
  
  # EE vs control
  lnRR_E2 <- log(EC_mean) - log(CC_mean)
  
  
  lnRRV_E2 <-  (EC_SD^2 / (EC_mean^2*EC_n)) + 
                            (CC_SD^2 / (CC_mean^2*CC_n))
  
  
  # Stress vs control
  lnRR_S2 <- log(CS_mean) - log(CC_mean)
  
  lnRRV_S2 <- (CS_SD^2 / (CS_mean^2*CS_n)) + 
                           (CC_SD^2 / (CC_mean^2*CC_n))
  
  # EE + stress vs control
  lnRR_ES2 <- log(ES_mean) - log(CC_mean)
  
  lnRRV_ES2 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                            (CC_SD^2 / (CC_mean^2*CC_n))
  
  # EE + stress vs stress (the effect of E in the presence of S)
  lnRR_E3 <- log(ES_mean) - log(CS_mean)
  
  lnRRV_E3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                           (CS_SD^2 / (CS_mean^2*CS_n))
  
  # EE + stress vs EE (the effect of S in the presence of E)
  lnRR_S3 <- log(ES_mean) - log(EC_mean)
  
  lnRRV_S3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                           (EC_SD^2 / (EC_mean^2*EC_n))
  
  effect2 <- tibble(
    lnRR_E2 = lnRR_E2,
    lnRRV_E2 = lnRRV_E2, 
    lnRR_S2 = lnRR_S2, 
    lnRRV_S2 = lnRRV_S2, 
    lnRR_ES2 =lnRR_ES2, 
    lnRRV_ES2 = lnRRV_ES2,
    lnRR_E3 =lnRR_E3, 
    lnRRV_E3 = lnRRV_E3,
    lnRR_S3 = lnRR_S3,
    lnRRV_S3 = lnRRV_S3
  )
  effect2
  }
  else {
    asin_trans <- function(percent) { asin(sqrt(percent/100)) }
    
    # transforming SD 
    ES_SD <- sqrt((ES_SD/100)^2/(4*(ES_mean/100)*(1-(ES_mean/100))))
    EC_SD <- sqrt((EC_SD/100)^2/(4*(EC_mean/100)*(1-(EC_mean/100))))
    CS_SD <- sqrt((CS_SD/100)^2/(4*(CS_mean/100)*(1-(CS_mean/100))))
    CC_SD <- sqrt((CC_SD/100)^2/(4*(CC_mean/100)*(1-(CC_mean/100))))
    
    # transformaing mean
    ES_mean <- asin_trans(ES_mean)
    EC_mean <- asin_trans(EC_mean)
    CS_mean <- asin_trans(CS_mean)
    CC_mean <- asin_trans(CC_mean)
    
    # EE vs control
    lnRR_E2 <- log(EC_mean) - log(CC_mean)
    
    
    lnRRV_E2 <- (EC_SD^2 / (EC_mean^2*EC_n)) + 
                              (CC_SD^2 / (CC_mean^2*CC_n))
    
    # Stress vs control
    lnRR_S2 <- log(CS_mean) - log(CC_mean)
    
    lnRRV_S2 <- (CS_SD^2 / (CS_mean^2*CS_n)) + 
                             (CC_SD^2 / (CC_mean^2*CC_n))
    
    # EE + stress vs control
    lnRR_ES2 <- log(ES_mean) - log(CC_mean)
    
    lnRRV_ES2 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                              (CC_SD^2 / (CC_mean^2*CC_n))
    
    # EE + stress vs stress (the effect of E in the presence of S)
    lnRR_E3 <-log(ES_mean) - log(CS_mean)
    
    lnRRV_E3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                             (CS_SD^2 / (CS_mean^2*CS_n))
    
    # EE + stress vs EE (the effect of S in the presence of E)
    lnRR_S3 <- log(ES_mean) - log(EC_mean)
    
    lnRRV_S3 <- (ES_SD^2 / (ES_mean^2*ES_n)) + 
                             (EC_SD^2 / (EC_mean^2*EC_n))
    
    
    effect2 <- tibble(
      lnRR_E2 = lnRR_E2,
      lnRRV_E2 = lnRRV_E2, 
      lnRR_S2 = lnRR_S2, 
      lnRRV_S2 = lnRRV_S2, 
      lnRR_ES2 =lnRR_ES2, 
      lnRRV_ES2 = lnRRV_ES2,
      lnRR_E3 =lnRR_E3, 
      lnRRV_E3 = lnRRV_E3,
      lnRR_S3 = lnRR_S3,
      lnRRV_S3 = lnRRV_S3
    )
    effect2
  }
  
}



#' @title make_VCV_matrix 
#' @description Function for generating simple covariance and correlation matrices based on a clustered variable
#' @param data Dataframe object containing effect sizes, their variance, unique IDs and clustering variable
#' @param matrix Sometimes clustering can get quite complicated. Here you can 'daisy' chain matrices for different levels of clustering to build a combined matrix. Just add in the matrix with an existing cluster and set up a new cluster.
#' @param V Name of the variable (as a string – e.g, "V1") containing effect size variances variances
#' @param m Mean of the control group that is shared. Only used when vcal does not equal "none".
#' @param sd Standard deviation of the control group that is shared. Only used when vcal does not equal "none".
#' @param n Sample size of the control group that is shared.Only used when vcal does not equal "none".
#' @param cluster Name of the variable (as a string – e.g, "V1") indicating which effects belong to the same cluster. Same value of 'cluster' are assumed to be nonindependent (correlated).
#' @param obs Name of the variable (as a string – e.g, "V1") containing individual IDs for each value in the V (Vector of variances). If this parameter is missing, label will be labelled with consecutive integers starting from 1.
#' @param rho Known or assumed correlation value among effect sizes sharing same 'cluster' value. Default value is 0.5.
#' @param type Optional logical parameter indicating whether a full variance-covariance matrix (default or "vcv") is needed or a correlation matrix ("cor") for the non-independent blocks of variance values.
#' @param vcal The calculation of the covariance. Defaults to "none" in which case rho is used. Otherwise, "ROM" (log response ratio) or "LOR" (log odds ratios) can be calculated based on large-sample approximations. 
#' @examples {
#' data(sparrows)
#' # Add sampling variance
#' sparrows$v <- 1 / (sparrows$SampleSize - 3)
#' # Fake grouping
#' sparrows$group1 <- rep(1:2, length.out = nrow(sparrows))
#' # Calculate V based on Place grouping, just for demonstration purposes. 
#' V <- make_VCV_matrix(data = sparrows, V = "v", cluster = "Place", type = "vcv", vcal = "none", rho = 0.5)
#' # Now we an add the second level of clustering
#' V <- make_VCV_matrix(data = sparrows, matrix = V, V = "v", cluster = "group1", type = "vcv", vcal = "none", rho = 0.5)
#' 
#' }
#' @export

make_VCV_matrix <- function(data, matrix = NULL, V, m, sd, n, cluster, obs, type=c("vcv", "cor"), vcal = c("none", "lnOR", "ROM"), rho=0.5){
  
  type <- match.arg(type, choices = type)
  vcal <- match.arg(vcal, choices = vcal)

  if (missing(data)) {
    stop("Must specify dataframe via 'data' argument.")
  }
  if (missing(V)) {
    stop("Must specify name of the variance variable via 'V' argument.")
  }
  if (missing(cluster)) {
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  }
  if (missing(obs)) {
    obs <- 1:length(V)   
  }
  if (missing(type)) {
    type <- "vcv" 
  }
  
  if (vcal != "none") {
    if(missing(m) | missing(n) | missing(sd)){
    stop("Must specify m, sd, and n arguments so that the covariance can be correctly calculated.")}
  }

  if(is.null(matrix)){
    new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) #make empty matrix of the same size as data length
    rownames(new_matrix) <- data[ ,obs]
    colnames(new_matrix) <- data[ ,obs]
    } else {
      new_matrix <- matrix
    }

  # find start and end coordinates for the subsets
  tmp <- duplicated(data[ ,cluster])
  shared_coord <- which(data[ ,cluster] %in% data[tmp, cluster]==TRUE)
  # matrix of combinations of coordinates for each experiment with shared control
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(utils::combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      
      if(vcal == "none"){
        p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      } else {
        p1_p2_cov <- vcalc(m, sd, n, p1, data, cov_type = vcal)
      }

      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]   #add the diagonal
  }
  
  if(type == "cor"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1   #add the diagonal of 1
  }
  
  rownames(new_matrix) <- 1:dim(new_matrix)[1]
  colnames(new_matrix) <- 1:dim(new_matrix)[1]

  return(new_matrix)
}

###########
# Function to calculate power (two-tail) for meta-analysis
power.ma_Shinichi <- function(mu, SE, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/SE)-pnorm(qnorm(1-alpha/2)+abs(mu)/SE)
} # or power.ma_Shinichi1 <- function(mu,SE){1 - pnorm(qnorm(1-0.05/2)-abs(mu)/SE) + pnorm(-qnorm(1-0.05/2)-abs(mu)/SE)}


# Function for power analysis for empirical data point
power.individual_Shinichi <- function(mu, se, alpha = 0.05) {
  2-pnorm(qnorm(1-alpha/2)-abs(mu)/se)-pnorm(qnorm(1-alpha/2)+abs(mu)/se)} # two-tailed power


# Function for Type S error for empirical data point
error_S <- function(mu, se, alpha = 0.05){
  #z <- qnorm(1 - alpha/2) # Z-score or quantile
  p.u <- 1 - pnorm(qnorm(1 - alpha/2) - abs(mu)/se) # upper-tail probability
  p.l <- pnorm(-qnorm(1 - alpha/2) - abs(mu)/se) # lower-tail probability
  power <- p.u + p.l # upper + lower
  errorS <- p.l/power # percentage of the opposite direction
  return(errorS)
} 

# Function for Type M error for empirical data point
error_M <- function(mu, se, alpha = 0.05, N = 10000) {
  est.random <- rnorm(n=N, mean = mu, sd = se)
  # est.random <- mu + se*rnorm(n=N, mean=0, sd=1)
  sig.index <- abs(est.random) > se*qnorm(1 - alpha/2)
  overestimate <- mean(abs(est.random)[sig.index])/abs(mu) # ratio is regardnesss of sign, so we need absolute value
  absolute_error <- overestimate*abs(mu) - abs(mu)
  relative_error <- absolute_error/(overestimate*abs(mu))
  return(abs(overestimate) %>% round(3))
}


# Get the results from the model
pred_dist_data <- function(mod_sim){
  # Get sigmas
  sigma2 <- mod_sim$sigma2
  
  # Get sampling error for mean
  se2_mu <- mod_sim$se^2
  
  # Get mean estimate
  mu <- mod_sim$b[1]
  
  # Names of random effect levels
  names <- c(attr(mod_sim$s.names, "names"), 'total')
  
  # Sigmas for prediction interval at each level
  cum_sigma2 <- c(sigma2 + se2_mu, sum(sigma2, se2_mu))
  
  # Create the dataframe
  data <- data.frame(group =  names, 
                     mean = mu,
                     sd = sqrt(cum_sigma2))
  return(data)
}

#### Calculate the proportion of effects beyond some threshold
prop_beyond <- function(data, threshold = 0.2, tail = c("above", "below")) {
  tail = match.arg(tail)
  
  # Get the mean and sd
  mean <- data$mean
  sd <- data$sd
  
  # Get the proportion beyond the threshold
  if (tail == "above") {
    proportion <- pnorm(threshold, mean, sd, lower.tail = FALSE)
  } else if (tail == "below") {
    proportion <- pnorm(threshold, mean, sd, lower.tail = TRUE)
  } else {
    stop("tail must be either 'above' or 'below'")
  }
  
  # Return the proportion
  return(round(proportion*100, 2))
}

#--------------------------t distribution--------------------------#
library(extraDistr)

#### Calculate the proportion of effects beyond some threshold
propT_beyond <- function(data, df, threshold = 0.2, tail = c("above", "below")) {
  tail = match.arg(tail)
  
  # Get the mean and sd
  m <- data$mean
  sd <- data$sd
  
  # Get the proportion beyond the threshold
  if (tail == "above") {
    proportion <- extraDistr::plst(threshold, df = df, mu = m, sigma = sd, lower.tail = FALSE)
  } else if (tail == "below") {
    proportion <- extraDistr::plst(threshold, df= df, mu = m, sigma = sd, lower.tail = TRUE)
  } else {
    stop("tail must be either 'above' or 'below'")
  }
  
  # Return the proportion
  return(round(proportion*100, 2))
}
