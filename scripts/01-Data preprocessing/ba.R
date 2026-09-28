#install.packages("remotes")
#remotes::install_github("rjaneUCF/MultiHazard")

#Loading libraries
library(MultiHazard)
library(VineCopula)

#watershed name
name= 'Mullica'

#Cora file
cora_path = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'_cora_centroid_detrend_ntr_ts.csv',sep="")
cora_path_decl = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'_cora_centroid_detrend_declust_ntr_ts.csv',sep="")

#Distances 
dist = 5

#Basin-average AORC totals
cora_labels = c("A", "B", "C", "D", "E")

# Putting all ba series for a given distance into a dataframe
ba_aorc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'_aorc_ba_ts',sep=""))

# Centroid of AORC points
ba_lat_lon = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'/Wahl Approach/'/name,'_aorc_centroid_km.csv', sep=""))

#years
years = 1979:2024

#Return period values at the top
rp_rain <- 10
rp_ntr <- 10

#Rainfall event duration
rain_dur = 24

# Helper functions (these need to be defined before parallel processing)
extract_co_max <- function(event_ids, data, window) {
  sapply(event_ids, function(id) {
    start_idx <- max(1, id - window/2)
    end_idx <- min(length(data), id + window/2)
    max(data[start_idx:end_idx], na.rm = TRUE)
  })
}


classify_event_types <- function(event_ids, names_data, start_idx, end_idx) {
  sapply(seq_along(event_ids), function(i) {
    event_range <- start_idx[i]:end_idx[i]
    ifelse(any(!is.na(names_data[event_range])), "tc", "non-tc")
  })
}


#Function to calculate correlations
calculate_correlations <- function(event_ids, type_ids, precip_data, wl_data) {
  tc_mask <- type_ids == "tc"
  non_tc_mask <- type_ids == "non-tc"
  
  cor_tc <- if(sum(tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[tc_mask], wl_data[tc_mask],
                              method="kendall")$estimate)
  } else NA
  
  cor_non_tc <- if(sum(non_tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[non_tc_mask], wl_data[non_tc_mask],
                              method="kendall")$estimate)
  } else NA
  
  list(tc = as.numeric(cor_tc), non_tc = as.numeric(cor_non_tc))
}

calculate_correlations_with_sig<- function(event_ids, type_ids, precip_data, wl_data) {
  tc_mask <- type_ids == "tc"
  non_tc_mask <- type_ids == "non-tc"
  
  cor_tc <- if(sum(tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[tc_mask], wl_data[tc_mask],
                              method="kendall")$estimate)
  } else NA
  
  cor_non_tc <- if(sum(non_tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[non_tc_mask], wl_data[non_tc_mask],
                              method="kendall")$estimate)
  } else NA
  
  cor_pval_tc <- if(sum(tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[tc_mask], wl_data[tc_mask],
                              method="kendall")$p.value)
  } else NA
  
  cor_pval_non_tc <- if(sum(non_tc_mask) > 1) {
    suppressWarnings(cor.test(precip_data[non_tc_mask], wl_data[non_tc_mask],
                              method="kendall")$p.value)
  } else NA
  
  
  list(tc = as.numeric(cor_tc),  pval_tc = as.numeric(cor_pval_tc), 
       non_tc = as.numeric(cor_non_tc), pval_non_tc = as.numeric(cor_pval_non_tc))
}

calculate_copula <- function(event_ids, type_ids, precip_data, wl_data) {
  tc_mask <- type_ids == "tc"
  non_tc_mask <- type_ids == "non-tc"
  
  cop_tc <- if(sum(tc_mask) > 1) {
    tryCatch(
      BiCopSelect(u1 = pobs(precip_data[tc_mask]),
                  u2 = pobs(wl_data[tc_mask]),
                  familyset = c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)),
      error = function(e) list(family = NA, par = NA, par2 = NA)
    )  
  } else {
    list(family = NA, par = NA, par2 = NA)
  }
  
  cop_non_tc <- if(sum(non_tc_mask) > 1) {
    tryCatch(
      BiCopSelect(u1 = pobs(precip_data[non_tc_mask]),
                  u2 = pobs(wl_data[non_tc_mask])),
      error = function(e) list(family = NA, par = NA, par2 = NA)
    )
  } else {
    list(family = NA, par = NA, par2 = NA)
  }
  
  list(tc = c(cop_tc$family, cop_tc$par, cop_tc$par2),
       non_tc = c(cop_non_tc$family, cop_non_tc$par, cop_non_tc$par2))
}


# Calculate non-exceedance probability for a specific water level
calc_gpd_nonexceed_prob <- function(x, threshold, scale, shape, lambda) {
  
  # Above threshold: use GPD
  y <- (x - threshold) / scale
  
  if(abs(shape) < 1e-10) {
    # Shape ≈ 0 (exponential case)
    prob_exceed_threshold <- exp(-y)
  } else {
    # General case
    prob_exceed_threshold <- (1 + shape * y)^(-1/shape)
  }
  
  prob <- 1 - prob_exceed_threshold
  
  return(prob)
}

Diag_Non_Con_AIC<-function(Data,Omit=NA,x_lab=NA,y_lim_min=0,y_lim_max=1){
  
  #Load Gamlss package
  if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("The 'gamlss.dist' package is required but not installed.")
  }
  GU <- get("GU", envir = asNamespace("gamlss.dist"))
  RG <- get("RG", envir = asNamespace("gamlss.dist"))
  
  # Load density and cumulative functions safely
  dRG <- get("dRG", envir = asNamespace("gamlss.dist"))
  pRG <- get("pRG", envir = asNamespace("gamlss.dist"))
  dGU <- get("dGU", envir = asNamespace("gamlss.dist"))
  pGU <- get("pGU", envir = asNamespace("gamlss.dist"))
  
  #Distributions to test
  Dist<-c("Gaus","Gum","Lapl","Logis","RGum")
  Dist.2 <- Dist
  Test<-1:5
  if(!is.na(Omit[1])){
    Test<-Test[-which(Dist %in% Omit)]
    Dist.2 <- Dist.2[-which(Dist %in% Omit)]
  }
  
  #Second omit vector for distributions where paramter estimation fails
  Omit.2 = rep(NA,3)
  
  par(mfrow=c(3,1))
  par(mar=c(4.2,4.2,1,1))
  
  #AIC result objects
  AIC.Gaus <- NA
  AIC.Gum <- NA
  AIC.Lapl <- NA
  AIC.Logis <- NA
  AIC.RGum <- NA
  
  #
  if(any(Test==2)){
    fit <- tryCatch(gamlss(Data  ~ 1, family=GU, trace=FALSE),
                    error = function(e) "error")
    if(fit[1] == "error") {
      Omit.2[1] = "Gum"
    } else {
      Omit.2[1] = ifelse(exp(fit$sigma.coefficients) < 0, "Gum", NA)
    }
  }
  
  if(any(Test==3)){
    fit <- tryCatch(fitdistr(Data, dlaplace, start=list(location=mean(Data),scale=sd(Data)/sqrt(2))),
                    error = function(e) "error")
    Omit.2[2] = ifelse(fit[1]=="error","Lapl",NA)
  }
  
  if(any(Test==5)){
    fit <- tryCatch(gamlss(Data ~ 1,family=RG, trace=FALSE),
                    error = function(e) "error")
    if(fit[1] == "error") {
      Omit.2[3] = "RGum"
    } else {
      Omit.2[3] = ifelse(exp(fit$sigma.coefficients) < 0, "RGum", NA)
    }
  }
  
  #Distributions to test
  if(any(is.na(Omit.2[1])==FALSE)){
    Test<-Test[-which(Dist.2 %in% Omit.2)]
  }
  
  #AIC
  if(any(Test==1)){
    fit<-fitdistr(Data, "normal")
    AIC.Gaus<-2*length(fit$estimate)-2*fit$loglik
  }
  
  if(any(Test==2)){
    fit <- gamlss(Data  ~ 1, family= GU, trace=FALSE)
    AIC.Gum <-fit$aic
  }
  
  if(any(Test==3)){
    fit <- fitdistr(Data, dlaplace, start=list(location=mean(Data), scale=sd(Data)/sqrt(2)))
    AIC.Lapl<-2*length(coef(fit))-2*logLik(fit)
  }
  
  
  if(any(Test==4)){
    fit<-fitdistr(Data,"logistic")
    AIC.Logis<-2*length(fit$estimate)-2*fit$loglik
  }
  
  if(any(Test==5)){
    fit <- gamlss(Data ~ 1,family=RG, trace=FALSE)
    AIC.RGum<-fit$aic
  }
  
  AIC<-data.frame(c("Gaus","Gum","Lapl","Logis","RGum"),c(AIC.Gaus,AIC.Gum,AIC.Lapl,AIC.Logis,AIC.RGum))
  colnames(AIC)<-c("Distribution","AIC")
  Best_fit<-AIC$Distribution[which(AIC$AIC==min(AIC$AIC,na.rm=T))]
  res<-list("AIC"=AIC,"Best_fit"=Best_fit)
  return(res)
}

Diag_Non_Con_Trunc_AIC<-function(Data,Omit=NA,x_lab="Data",y_lim_min=0,y_lim_max=1){
  
  # Validate Omit parameter
  valid_distributions <- c("BS", "Exp", "Gam(2)", "Gam(3)", "GamMix(2)", "GamMix(3)", "LNorm", "TNorm", "Twe", "Weib")
  if (!is.na(Omit[1])) {
    if (!all(Omit %in% valid_distributions)) {
      invalid_omit <- Omit[!Omit %in% valid_distributions]
      stop("Invalid distribution names in Omit.")
    }
    if (length(Omit) >= length(valid_distributions)) {
      stop("Cannot omit all distributions. At least one distribution must be tested.")
    }
  }
  
  #Checks gamlss package is installed
  if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
    stop("The 'gamlss.dist' package is required but not installed.")
  }
  GG <- gamlss.dist::GG #get("GG", envir = asNamespace("gamlss.dist"))
  GA <- get("GA", envir = asNamespace("gamlss.dist"))
  
  # Load density and cumulative functions safely
  dGG <- get("dGG", envir = asNamespace("gamlss.dist"))
  pGG <- get("pGG", envir = asNamespace("gamlss.dist"))
  dGA <- get("dGA", envir = asNamespace("gamlss.dist"))
  pGA <- get("pGA", envir = asNamespace("gamlss.dist"))
  
  #Distributions to test
  Dist<-c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")
  Test<-1:10
  if(is.na(Omit[1])==F){
    Test<-Test[-which(Dist %in% Omit)]
  }
  
  #AIC result objects
  AIC.BS<-NA
  AIC.Exp<-NA
  AIC.Gam2<-NA
  AIC.Gam3<-NA
  AIC.GamMix2<-NA
  AIC.GamMix3<-NA
  AIC.logNormal<-NA
  AIC.TNormal<-NA
  AIC.Tweedie<-NA
  AIC.Weib<-NA
  
  #AIC
  if(any(Test==1)){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = Data)
    fit.BS <- vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    AIC.BS<-2*length(coef(fit.BS))-2*logLik(fit.BS)
  }
  if(any(Test==2)){
    fit.Exp<-fitdistr(Data,"exponential")
    AIC.Exp<-2*length(fit.Exp$estimate)-2*fit.Exp$loglik
  }
  if(any(Test==3)){
    fit.Gam2<-fitdistr(Data, "gamma")
    AIC.Gam2<-2*length(fit.Gam2$estimate)-2*fit.Gam2$loglik
  }
  data.gamlss <- data.frame(X=Data)
  if(any(Test==4)){
    ### 3-parameter gamma dist.
    for(i in 1:100){
      fit.Gamma3 <- tryCatch(gamlss(X~1, data=data.gamlss, family=GG, trace=FALSE),
                             error = function(e) "error")
      if( is.character(fit.Gamma3) ) next
      if( !is.character(fit.Gamma3) ) break
    }
    if( is.character(fit.Gamma3) ){
      #AIC.Gamma3 <- -9999
      Test <- Test[-which(Test==4)]
    }else{
      AIC.Gam3 <- fit.Gamma3$aic
    }
  }
  if(any(Test==5)){
    ### 2 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX2_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=2, trace=FALSE),
                                 error = function(e) "error")
      if( is.character(fit.GamMIX2_GA) ) next
      if( !is.character(fit.GamMIX2_GA) ) break
    }
    if( is.character(fit.GamMIX2_GA) ){
      #AIC.GamMIX2_GA <- -9999
      Test <- Test[-which(Test==5)]
    }else{
      AIC.GamMix2 <- fit.GamMIX2_GA$aic
    }
  }
  if(any(Test==6)){
    ### 3 mixture-gamma dist.
    for(i in 1:100){
      fit.GamMIX3_GA <- tryCatch(gamlssMX(X~1, data=data.gamlss, family=GA, K=3, trace=FALSE),
                                 error = function(e) "error")
      if( is.character(fit.GamMIX3_GA) ) next
      if( !is.character(fit.GamMIX3_GA) ) break
    }
    if( is.character(fit.GamMIX3_GA) ){
      #AIC.GamMIX3_GA <- -9999
      Test <- Test[-which(Test==6)]
    }else{
      AIC.GamMix3 <- fit.GamMIX3_GA$aic
    }
  }
  #fit<-fitdist(Data, "invgauss", start = list(mean = 5, shape = 1))
  #AIC.InverseNormal<-2*length(fit$estimate)-2*fit$loglik
  if(any(Test==7)){
    fit.LNorm<-fitdistr(Data,"lognormal")
    AIC.logNormal<-2*length(fit.LNorm$estimate)-2*fit.LNorm$loglik
  }
  if(any(Test==8)){
    fit.TNorm <- fitdistr(Data, "normal")
    AIC.TNormal <- 2 * length(fit.TNorm$estimate) - 2 * fit.TNorm$loglik
  }
  if(any(Test==9)){
    capture.output(
      fit.Twe <- tweedie.profile(Data ~ 1,
                                 p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE),
      type = "output"
    )
    AIC.Tweedie<-2*3-2*fit.Twe$L.max
  }
  if(any(Test==10)){
    fit.Weib<-fitdistr(Data,"weibull")
    AIC.Weib<-2*length(fit.Weib$estimate)-2*fit.Weib$loglik
  }
  
  AIC<-data.frame(c("BS","Exp","Gam(2)","Gam(3)","GamMix(2)","GamMix(3)","LNorm","TNorm","Twe","Weib")[Test],c(AIC.BS,AIC.Exp,AIC.Gam2,AIC.Gam3,AIC.GamMix2,AIC.GamMix3,AIC.logNormal,AIC.TNormal,AIC.Tweedie,AIC.Weib)[Test])
  colnames(AIC)<-c("Distribution","AIC")
  Best_fit<-AIC$Distribution[which(AIC$AIC==min(AIC$AIC))]
  res<-list("AIC"=AIC, "Best_fit"=Best_fit)
  return(res)
}



joint_rp_func <- function(ba_aorc, ba_lat_lon, cora_labels, years, cora_path, cora_path_decl, name, rp_rain, rp_ntr, rain_dur){
  tryCatch({
    
    # AORC data (removing) last two years' worth of data as in original code
    all_aorc_time <- seq(from = as.POSIXct("1979-02-01 00:00:00", tz = "UTC"), to   = as.POSIXct("2022-12-31 23:00:00", tz = "UTC"),  by   = "hour")
    
    # Create dataframe
    aorc_precip_df <- data.frame(all_aorc_time, ba_aorc[1:(length(ba_aorc) - 366*24 - 365*24)])
    colnames(aorc_precip_df) <- c("date", "precip")
    
    #Find TC anc non-TC events rainfall volume
    aorc_precip_df <- HURDAT(aorc_precip_df,
                             lat.loc = ba_lat_lon$lat[i],
                             lon.loc = -ba_lat_lon$lon[i],
                             rad = 350)
    
    # Decluster events
    aorc_precip_df_decl <- Decluster_S_SW(Data = aorc_precip_df[,1:2],
                                          Window_Width_Sum = rain_dur,
                                          Window_Width = 5*24)
    print(summary(aorc_precip_df_decl))
    # Event identification
    ids <- order(aorc_precip_df_decl$Declustered, decreasing = TRUE)[1:(5*44)]
    types <- classify_event_types(ids, aorc_precip_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(aorc_precip_df), ids + 24))
    rain_event_data = list(ids = ids, types = types)
    
    #Fit GPD
    gpd_aorc = GPD_Fit(Data=aorc_precip_df_decl$Declustered[rain_event_data$ids], Data_Full=aorc_precip_df_decl$Declustered, Thres=min(aorc_precip_df_decl$Declustered[rain_event_data$ids]), Method = "Solari")
    gpd_aorc$Rate = length(rain_event_data$ids)/(nrow(aorc_precip_df)/(24*365.25))
    
    #Extract ten year event
    if (gpd_aorc$xi != 0) {
      rain_10yr <- gpd_aorc$Threshold + (gpd_aorc$sigma / gpd_aorc$xi) * ((rp_rain*gpd_aorc$Rate)^gpd_aorc$xi - 1)
    } else {
      rain_10yr <- gpd_aorc$Threshold + gpd_aorc$sigma * log(rp_rain*gpd_aorc$Rate)
    }
    
    #Early exit
    # Process each CORA location
    results_by_cora <- list()
    
    l = which(c("A","B","C","D","E") == cora_labels)
    print(l)
    #Read water level time series
    cora_ntr_df = read.csv(cora_path)[,-1][,c(1,l+1)]
    colnames(cora_ntr_df) <- c("date", "ntr")
    cora_ntr_df$date = as.POSIXct(cora_ntr_df$date)
    
    #Read in declustered water levels
    cora_ntr_decl_df = read.csv(cora_path_decl)[,-1][,c(1,l+1)]
    colnames(cora_ntr_decl_df) <- c("date", "ntr")
    cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)
    
    cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                               lat.loc = aorc_subset_lat_lon$lat[i],
                               lon.loc = -aorc_subset_lat_lon$lon[i],
                               rad = 350)
    
    # Event identification
    print(summary(cora_ntr_decl_df))
    ids <- order(cora_ntr_decl_df$ntr, decreasing = TRUE)[1:(5*44)]
    types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 24))
    ntr_event_data = list(ids = ids, types = types)
    
    evt = ntr_event_data
    ntr_ext = evt$ids
    
    #Fit GPD
    gpd_ntr = GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext]), Method="Solari")
    gpd_ntr$Rate = length(cora_ntr_decl_df$ntr[ntr_ext])/(nrow(cora_ntr_df)/(24*365.25))
    
    #Extract ten year event
    if (gpd_ntr$xi != 0) {
      ntr_10yr <- gpd_ntr$Threshold + (gpd_ntr$sigma / gpd_ntr$xi) * ((rp_ntr*gpd_ntr$Rate)^gpd_ntr$xi - 1)
    } else {
      ntr_10yr <- gpd_ntr$Threshold + gpd_ntr$sigma * log(rp_ntr*gpd_ntr$Rate)
    }
    
    v_ntr <- extract_co_max(event_ids=ntr_ext, data=aorc_precip_df_decl$Totals, window=3*24) 
    
    #Correlation without stratification by generating mechanism
    cor = cor.test(v_ntr,cora_ntr_df$ntr[ntr_ext], method="kendall")
    
    cor_ntr <- cor$estimate
    cor_pval_ntr <- cor$p.value
    
    #Characteristics of sample con. on ntr
    correlations <- calculate_correlations_with_sig(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])
    
    cor_ntr_tc <- as.numeric(correlations$tc)
    cor_ntr_non_tc <- as.numeric(correlations$non_tc)
    
    cor_pval_ntr_tc <- as.numeric(correlations$pval_tc)
    cor_pval_ntr_non_tc <- as.numeric(correlations$pval_non_tc)
    
    copula = calculate_copula(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])
    
    cop_ntr_tc_family <- as.numeric(copula$tc[1])
    cop_ntr_non_tc_family <- as.numeric(copula$non_tc[1])
    
    cop_ntr_tc_par1 <- as.numeric(copula$tc[2])
    cop_ntr_non_tc_par1 <- as.numeric(copula$non_tc[2])
    
    cop_ntr_tc_par2 <- as.numeric(copula$tc[3])
    cop_ntr_non_tc_par2 <- as.numeric(copula$non_tc[3])
    
    tc_indices <- evt$types == "tc"
    non_tc_indices <- evt$types == "non-tc"
    
    n_ntr_tc = sum(tc_indices)
    n_ntr_non_tc = sum(non_tc_indices)
    
    cop_ntr_tc =BiCopSelect(u1 = pobs(v_ntr[tc_indices]),
                            u2 = pobs(cora_ntr_df$ntr[ntr_ext][tc_indices]),
                            familyset = c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36))
    
    cop_ntr_non_tc =BiCopSelect(u1 = pobs(v_ntr[non_tc_indices]),
                                u2 = pobs(cora_ntr_df$ntr[ntr_ext][non_tc_indices]))
    
    #Fit GPDs
    gpd_ntr_tc =       GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext][tc_indices], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext][tc_indices]), Method = "Solari")
    gpd_ntr_non_tc =   GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext][non_tc_indices], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext][non_tc_indices]), Method = "Solari")
    
    gpd_ntr_tc$Rate =  n_ntr_tc/(nrow(aorc_precip_df)/(24*365.25))
    gpd_ntr_non_tc$Rate =  n_ntr_non_tc/(nrow(aorc_precip_df)/(24*365.25))
    
    #Find u of ntr_10yr
    ntr_10yr_tc_u = calc_gpd_nonexceed_prob(ntr_10yr, gpd_ntr_tc$Threshold, gpd_ntr_tc$sigma, gpd_ntr_tc$xi, lambda=gpd_ntr_tc$Rate)
    ntr_10yr_non_tc_u =  calc_gpd_nonexceed_prob(ntr_10yr, gpd_ntr_non_tc$Threshold, gpd_ntr_non_tc$sigma, gpd_ntr_non_tc$xi, lambda=gpd_ntr_non_tc$Rate)
    
    v_ntr_tc = v_ntr[tc_indices] + runif(length(v_ntr[tc_indices]),0.0001,0.001)
    v_ntr_non_tc = v_ntr[non_tc_indices] + runif(length(v_ntr[non_tc_indices]),0.0001,0.001)
    
    #Find best fitting distribution for the conditiionec variable
    non_con_dist_r_tc = Diag_Non_Con_Trunc_AIC(v_ntr_tc, Omit="Weib")$Best_fit
    non_con_dist_r_non_tc = Diag_Non_Con_Trunc_AIC(v_ntr_non_tc, Omit="Weib")$Best_fit
    
    #Marginal distributions for rainfall
    #Tc
    if(non_con_dist_r_tc=="BS"){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = v_ntr_tc)
      non_con_dist_fit_r_tc<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
      r_10yr_con_ntr_tc_u<-pbisa(rain_10yr, as.numeric(Coef(non_con_dist_fit_r_tc)[1]),as.numeric(Coef(non_con_dist_fit_r_tc)[2]))
    }
    if(non_con_dist_r_tc=="Exp"){
      non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc,"exponential")
      r_10yr_con_ntr_tc_u<-pexp(rain_10yr, as.numeric(non_con_dist_fit_r_tc$estimate[1]))
    }
    if(non_con_dist_r_tc=="Gam(2)"){
      non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc, "gamma")
      r_10yr_con_ntr_tc_u<-pgamma(rain_10yr, shape = as.numeric(non_con_dist_fit_r_tc$estimate[1]), rate = as.numeric(non_con_dist_fit_r_tc$estimate[2]))
    }
    if(non_con_dist_r_tc=="Gam(3)"){
      data.gamlss = data.frame(X=v_ntr_tc)
      non_con_dist_fit_r_tc <-  tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                         error = function(e) "error")
      r_10yr_con_ntr_tc_u<-pGG(rain_10yr, mu=exp(non_con_dist_fit_r_tc$mu.coefficients), sigma=exp(non_con_dist_fit_r_tc$sigma.coefficients), nu=non_con_dist_fit_r_tc$nu.coefficients)
    }
    if(non_con_dist_r_tc=="InvG"){
      non_con_dist_fit_r_tc <- fitdist(v_ntr_tc, "invgauss", start = list(mean = 5, shape = 1))
      r_10yr_con_ntr_tc_u<-pinvgauss(rain_10yr, as.numeric(non_con_dist_fit_r_tc$estimate[1]), as.numeric(non_con_dist_fit_r_tc$estimate[2]))
    }
    if(non_con_dist_r_tc=="LogN"){
      non_con_dist_fit_r_tc <- fitdistr(v_ntr_tc,"lognormal")
      r_10yr_con_ntr_tc_u<-plnorm(rain_10yr, meanlog = as.numeric(non_con_dist_fit_r_tc$estimate[1]), sdlog = as.numeric(non_con_dist_fit_r_tc$estimate[2]))
    }
    if(non_con_dist_r_tc=="TNorm"){
      non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc,"normal")
      r_10yr_con_ntr_tc_u<-ptruncnorm(rain_10yr,a=min(v_ntr_tc),as.numeric(non_con_dist_fit_r_tc$estimate[1]),as.numeric(non_con_dist_fit_r_tc$estimate[2]))
    }
    if(non_con_dist_r_tc=="Twe"){
      non_con_dist_fit_r_tc <-tweedie.profile(v_ntr_tc ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
      r_10yr_con_ntr_tc_u<-ptweedie(rain_10yr, power=non_con_dist_fit_r_tc$p.max, mu=mean(v_ntr_tc), phi=non_con_dist_fit_r_tc$phi.max)
    }
    if(non_con_dist_r_tc=="Weib"){
      non_con_dist_fit_r_tc <- fitdistr(v_ntr_tc, "weibull")
      r_10yr_con_ntr_tc_u<-pweibull(rain_10yr, as.numeric(non_con_dist_fit_r_tc$estimate[1]), as.numeric(non_con_dist_fit_r_tc$estimate[2]))
    }
    
    #non-Tc
    if(non_con_dist_r_non_tc=="BS"){
      bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
      bdata2 <- transform(bdata2, y = v_ntr_non_tc)
      non_con_dist_fit_r_non_tc<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
      r_10yr_con_ntr_non_tc_u<-pbisa(rain_10yr, as.numeric(Coef(non_con_dist_fit_r_non_tc)[1]),as.numeric(Coef(non_con_dist_fit_r_non_tc)[2]))
    }
    if(non_con_dist_r_non_tc=="Exp"){
      non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc,"exponential")
      r_10yr_con_ntr_non_tc_u<-pexp(rain_10yr, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]))
    }
    if(non_con_dist_r_non_tc=="Gam(2)"){
      non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc, "gamma")
      r_10yr_con_ntr_non_tc_u<-pgamma(rain_10yr, shape = as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), rate = as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
    }
    if(non_con_dist_r_non_tc=="Gam(3)"){
      data.gamlss = data.frame(X=v_ntr_non_tc)
      non_con_dist_fit_r_non_tc <-  tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                             error = function(e) "error")
      r_10yr_con_ntr_non_tc_u<-pGG(rain_10yr, mu=exp(non_con_dist_fit_r_non_tc$mu.coefficients), sigma=exp(non_con_dist_fit_r_non_tc$sigma.coefficients), nu=non_con_dist_fit_r_non_tc$nu.coefficients)
    }
    if(non_con_dist_r_non_tc=="InvG"){
      non_con_dist_fit_r_non_tc <- fitdist(v_ntr_non_tc, "invgauss", start = list(mean = 5, shape = 1))
      r_10yr_con_ntr_non_tc_u<-pinvgauss(rain_10yr, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
    }
    if(non_con_dist_r_non_tc=="LogN"){
      non_con_dist_fit_r_non_tc <- fitdistr(v_ntr_non_tc,"lognormal")
      r_10yr_con_ntr_non_tc_u<-plnorm(rain_10yr, meanlog = as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), sdlog = as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
    }
    if(non_con_dist_r_non_tc=="TNorm"){
      non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc,"normal")
      r_10yr_con_ntr_non_tc_u<-ptruncnorm(rain_10yr,a=min(v_ntr_non_tc),as.numeric(non_con_dist_fit_r_non_tc$estimate[1]),as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
    }
    if(non_con_dist_r_non_tc=="Twe"){
      non_con_dist_fit_r_non_tc <-tweedie.profile(v_ntr_non_tc ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
      r_10yr_con_ntr_non_tc_u<-ptweedie(rain_10yr, power=non_con_dist_fit_r_non_tc$p.max, mu=mean(v_ntr_non_tc), phi=non_con_dist_fit_r_non_tc$phi.max)
    }
    if(non_con_dist_r_non_tc=="Weib"){
      non_con_dist_fit_r_non_tc <- fitdistr(v_ntr_non_tc, "weibull")
      r_10yr_con_ntr_non_tc_u<-pweibull(rain_10yr, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
    }
    
    #Find non-exceedance probability
    #Obtaining sample conditioned on rainfall
    evt <- rain_event_data
    ntr_vals <- extract_co_max(evt$ids, cora_ntr_df$ntr, window=3*24)
    
    cat("ntr_vals:", ntr_vals, "\n")
    
    #Correlation without stratification by generating mechanism
    cor = cor.test(aorc_precip_df_decl$Declustered[evt$ids],ntr_vals, method="kendall")
    cor_r <- cor$estimate
    cor_pval_r <- cor$p.value
    
    #Characteristics of sample con. on r
    correlations <- calculate_correlations_with_sig(evt$ids, evt$types,
                                                    aorc_precip_df_decl$Declustered[evt$ids],
                                                    ntr_vals)
    
    cor_r_tc <- as.numeric(correlations$tc)
    cor_r_non_tc <- as.numeric(correlations$non_tc)
    
    cor_pval_r_tc <- as.numeric(correlations$pval_tc)
    cor_pval_r_non_tc <- as.numeric(correlations$pval_non_tc)
    
    copula = calculate_copula(evt$ids, evt$types,
                              aorc_precip_df_decl$Declustered[evt$ids],
                              ntr_vals)
    
    cop_r_tc_family <- as.numeric(copula$tc[1])
    cop_r_non_tc_family <- as.numeric(copula$non_tc[1])
    
    cop_r_tc_par1 <- as.numeric(copula$tc[2])
    cop_r_non_tc_par1 <- as.numeric(copula$non_tc[2])
    
    cop_r_tc_par2 <- as.numeric(copula$tc[3])
    cop_r_non_tc_par2 <- as.numeric(copula$non_tc[3])
    
    tc_r_indices <- evt$types == "tc"
    non_tc_r_indices <- evt$types == "non-tc"
    
    n_r_tc = sum(tc_r_indices)
    n_r_non_tc = sum(non_tc_r_indices)
    
    
    cop_r_tc <- BiCopSelect(u1 = pobs(aorc_precip_df_decl$Declustered[evt$ids][tc_r_indices]),
                            u2 = pobs(ntr_vals[tc_r_indices]),
                            familyset = c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36))
    
    cop_r_non_tc  <- BiCopSelect(u1 = pobs(aorc_precip_df_decl$Declustered[evt$ids][non_tc_r_indices]),
                                 u2 = pobs(ntr_vals[non_tc_r_indices]))
    
    
    #Fit GPDs
    gpd_rain_tc = GPD_Fit(Data=aorc_precip_df_decl$Declustered[evt$ids][tc_r_indices], Data_Full=aorc_precip_df_decl$Declustered[evt$ids][tc_r_indices], Thres=min(aorc_precip_df_decl$Declustered[evt$ids][tc_r_indices]), Method="Solari")
    gpd_rain_non_tc = GPD_Fit(Data=aorc_precip_df_decl$Declustered[evt$ids][non_tc_r_indices], Data_Full=aorc_precip_df_decl$Declustered[evt$ids][non_tc_r_indices], Thres=min(aorc_precip_df_decl$Declustered[evt$ids][non_tc_r_indices]), Method="Solari")
    
    #Correcting the rate
    gpd_rain_tc$Rate = sum(tc_r_indices)/(nrow(aorc_precip_df)/(24*365.25))
    gpd_rain_non_tc$Rate =  sum(non_tc_r_indices)/(nrow(aorc_precip_df)/(24*365.25))
    
    
    #Find u of ntr_10yr
    r_10yr_tc_u = calc_gpd_nonexceed_prob(rain_10yr, gpd_rain_tc$Threshold, gpd_rain_tc$sigma, gpd_rain_tc$xi, lambda=gpd_rain_tc$Rate)
    r_10yr_non_tc_u =  calc_gpd_nonexceed_prob(rain_10yr, gpd_rain_non_tc$Threshold, gpd_rain_non_tc$sigma, gpd_rain_non_tc$xi, lambda=gpd_rain_non_tc$Rate)
    
    non_con_dist_ntr_tc = Diag_Non_Con_AIC(ntr_vals[tc_r_indices],Omit="Lapl")$Best_fit
    non_con_dist_ntr_non_tc = Diag_Non_Con_AIC(ntr_vals[non_tc_r_indices],Omit="Lapl")$Best_fit
    
    #Marginal distributions for ntr
    #Tc
    if(non_con_dist_ntr_tc == "Gum"){
      non_con_dist_fit_ntr_tc <- gamlss(ntr_vals[tc_r_indices]  ~ 1, family= GU)
      ntr_10yr_con_r_tc_u<-pGU(ntr_10yr,as.numeric(non_con_dist_fit_ntr_tc$mu.coefficients),exp(as.numeric(non_con_dist_fit_ntr_tc$sigma.coefficients)))
    }
    if(non_con_dist_ntr_tc=="RGum"){
      non_con_dist_fit_ntr_tc <- gamlss(ntr_vals[tc_r_indices] ~ 1,family=RG)
      ntr_10yr_con_r_tc_u<-pRG(ntr_10yr,non_con_dist_fit_ntr_tc$mu.coefficients,exp(non_con_dist_fit_ntr_tc$sigma.coefficients))
    }
    if(non_con_dist_ntr_tc=="Gaus"){
      non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices],"normal")
      ntr_10yr_con_r_tc_u<-pnorm(ntr_10yr, as.numeric(non_con_dist_fit_ntr_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
    }
    if(non_con_dist_ntr_tc=="Lapl"){
      non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices], dlaplace, start=list(location=mean(ntr_vals[tc_r_indices]), scale=sd(ntr_vals[tc_r_indices])/sqrt(2)))
      ntr_10yr_con_r_tc_u <- plaplace(ntr_10yr,as.numeric(non_con_dist_fit_ntr_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
    }
    if(non_con_dist_ntr_tc=="Logis"){
      non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices],"logistic")
      ntr_10yr_con_r_tc_u<-plogis(ntr_10yr,as.numeric(non_con_dist_fit_ntr_tc$estimate[1]),as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
    }
    
    
    #Non-Tc
    if(non_con_dist_ntr_non_tc == "Gum"){
      non_con_dist_fit_ntr_non_tc <- gamlss(ntr_vals[non_tc_r_indices]  ~ 1, family= GU)
      ntr_10yr_con_r_non_tc_u<-pGU(ntr_10yr,as.numeric(non_con_dist_fit_ntr_non_tc$mu.coefficients),exp(as.numeric(non_con_dist_fit_ntr_non_tc$sigma.coefficients)))
    }
    if(non_con_dist_ntr_non_tc=="RGum"){
      non_con_dist_fit_ntr_non_tc <- gamlss(ntr_vals[non_tc_r_indices] ~ 1,family=RG)
      ntr_10yr_con_r_non_tc_u<-pRG(ntr_10yr,non_con_dist_fit_ntr_non_tc$mu.coefficients,exp(non_con_dist_fit_ntr_non_tc$sigma.coefficients))
    }
    if(non_con_dist_ntr_non_tc=="Gaus"){
      non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices],"normal")
      ntr_10yr_con_r_non_tc_u<-pnorm(ntr_10yr, as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
    }
    if(non_con_dist_ntr_non_tc=="Lapl"){
      non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices], dlaplace, start=list(location=mean(ntr_vals[non_tc_r_indices]), scale=sd(ntr_vals[non_tc_r_indices])/sqrt(2)))
      ntr_10yr_con_r_non_tc_u <- plaplace(ntr_10yr,as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
    }
    if(non_con_dist_ntr_non_tc=="Logis"){
      non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices],"logistic")
      ntr_10yr_con_r_non_tc_u<-plogis(ntr_10yr,as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]),as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
    }
    
    #Evaluate the copulas CDF
    UU_r_tc<-BiCopCDF(r_10yr_tc_u, ntr_10yr_con_r_tc_u, cop_r_tc)
    UU_r_non_tc<-BiCopCDF(r_10yr_non_tc_u, ntr_10yr_con_r_non_tc_u, cop_r_non_tc)
    UU_ntr_tc<-BiCopCDF(r_10yr_con_ntr_tc_u, ntr_10yr_tc_u, cop_ntr_tc)
    UU_ntr_non_tc<-BiCopCDF(r_10yr_con_ntr_non_tc_u, ntr_10yr_non_tc_u, cop_ntr_non_tc)
    
    
    #Interarrival time of events in sample conditioning on var 2
    EL_r_tc<-1/(n_r_tc/44)
    EL_r_non_tc<-1/(n_r_non_tc/44)
    EL_ntr_tc<-1/(n_ntr_tc/44)
    EL_ntr_non_tc<-1/(n_ntr_non_tc/44)
    
    #Function for evaluating the AEP at each point in the event space
    AEP_con_r_tc<-(1-r_10yr_tc_u-ntr_10yr_con_r_tc_u+UU_r_tc) / EL_r_tc
    AEP_con_r_non_tc <- (1-r_10yr_non_tc_u-ntr_10yr_con_r_non_tc_u+UU_r_non_tc) / EL_r_non_tc
    AEP_con_ntr_tc<-(1-r_10yr_con_ntr_tc_u-ntr_10yr_tc_u+UU_ntr_tc) / EL_ntr_tc
    AEP_con_ntr_non_tc <- (1-r_10yr_con_ntr_non_tc_u-ntr_10yr_non_tc_u+UU_ntr_non_tc) / EL_ntr_non_tc
    
    #Put the two AEPs into a dataframe
    AEP_tc = data.frame(AEP_con_r_tc,AEP_con_ntr_tc)
    AEP_non_tc = data.frame(AEP_con_r_non_tc,AEP_con_ntr_non_tc)
    
    #Select the maximum of the AEPs from the two samples at each grid point
    AEP_tc_max = apply(AEP_tc,1,function(x) max(x, na.rm = TRUE))
    AEP_non_tc_max = apply(AEP_non_tc,1,function(x) max(x, na.rm = TRUE))
    
    #Populations are independent so overall non-exceedance probability is product of individual non-exceedence probabilities
    ANEP = (1 - AEP_tc_max) * (1 - AEP_non_tc_max)
    
    #Convert non-excedence probabilities to return periods
    rp_joint = 1 / ( 1 - ANEP )
    
    #Exceedence probs combined with rate
    p_exceed_r_tc = 1-r_10yr_tc_u-ntr_10yr_con_r_tc_u+UU_r_tc
    p_exceed_ntr_tc = 1-r_10yr_con_ntr_tc_u-ntr_10yr_tc_u+UU_ntr_tc
    p_exceed_r_non_tc = 1-r_10yr_non_tc_u-ntr_10yr_con_r_non_tc_u+UU_r_non_tc
    p_exceed_ntr_non_tc = 1-r_10yr_con_ntr_non_tc_u-ntr_10yr_non_tc_u+UU_ntr_non_tc
    
    joint_annual_exceed_tc = max((n_r_tc/44)*p_exceed_r_tc,(n_ntr_tc/44)*p_exceed_ntr_tc)
    joint_annual_exceed_non_tc = max((n_r_non_tc/44)*p_exceed_r_non_tc,(n_ntr_non_tc/44)*p_exceed_ntr_non_tc)
    
    #Joint return period Poisson survival function
    rp_joint_poisson = 1 / (1 - exp(-(joint_annual_exceed_tc + joint_annual_exceed_non_tc)))
    
    # Store results
    res <- data.frame(
      cora = l,
      cor_r, cor_ntr,
      cor_r_tc = cor_r_tc, cor_r_non_tc = cor_r_non_tc,
      cor_ntr_tc = cor_ntr_tc, cor_ntr_non_tc = cor_ntr_non_tc,
      cor_pval_r = cor_pval_r, cor_pval_ntr = cor_pval_ntr,
      cor_pval_r_tc = cor_pval_r_tc, cor_pval_r_non_tc = cor_pval_r_non_tc,
      cor_pval_ntr_tc = cor_pval_ntr_tc, cor_pval_ntr_non_tc = cor_pval_ntr_non_tc,
      cop_r_tc_family = cop_r_tc_family, cop_r_non_tc_family = cop_r_non_tc_family,
      cop_ntr_tc_family  = cop_ntr_tc_family, cop_ntr_non_tc_family = cop_ntr_non_tc_family,
      cop_r_tc_par1 = cop_r_tc_par1 , cop_r_non_tc_par1 = cop_r_non_tc_par1,
      cop_ntr_tc_par1  = cop_ntr_tc_par1,cop_ntr_non_tc_par1 = cop_ntr_non_tc_par1,
      cop_r_tc_par2 = cop_r_tc_par2, cop_r_non_tc_par2 = cop_r_non_tc_par2,
      cop_ntr_tc_par2  = cop_ntr_tc_par2, cop_ntr_non_tc_par2 = cop_ntr_non_tc_par2,
      n_r_tc = n_r_tc, n_r_non_tc = n_r_non_tc,
      n_ntr_tc = n_ntr_tc, n_ntr_non_tc  = n_ntr_non_tc,
      non_con_dist_r_tc = non_con_dist_r_tc, non_con_dist_r_non_tc = non_con_dist_r_non_tc,
      non_con_dist_ntr_tc  = non_con_dist_ntr_tc, non_con_dist_ntr_non_tc = non_con_dist_ntr_non_tc,
      AEP_con_r_tc = AEP_con_r_tc, AEP_con_r_non_tc = AEP_con_r_non_tc,
      AEP_con_ntr_tc = AEP_con_ntr_tc, AEP_con_ntr_non_tc = AEP_con_ntr_non_tc,
      ANEP = ANEP,
      rp_joint,
      p_exceed_r_tc, p_exceed_ntr_tc,
      p_exceed_r_non_tc, p_exceed_ntr_non_tc,
      joint_annual_exceed_tc, joint_annual_exceed_non_tc,
      rp_joint_poisson,
      rain_10yr = rain_10yr, ntr_10yr = ntr_10yr
    )
    print(res)
    res
    
  }, error = function(e) {
    message(paste("Error processing location", e$message))
    traceback()  # Add this
    return(NULL)
  })
}

result_df = joint_rp_func(ba_aorc = ba_aorc, ba_lat_lon = ba_lat_lon,
                          cora_labels =cora_labels, 
                          years=years, cora_path=cora_path, cora_path_decl=cora_path_decl, 
                          name=name, 
                          rp_rain=rp_rain, rp_ntr=rp_ntr, 
                          rain_dur=rain_dur)

# Write to CSV
write.csv(result_df, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'/Wahl Approach/',name,'_anv_spatial_',rain_dur,'_hr_acc_time_10_yr_r_ntr_cora_',cora_labels[1],'_dist_',dist,'.csv',sep=""))
