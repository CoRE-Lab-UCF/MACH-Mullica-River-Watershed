# Load required packages
library(MultiHazard)
library(VineCopula)
library(MASS)
library(tweedie)
library(gamlss)
library(gamlss.mx)
library(stats4)
library(copula)
library(VGAM)
library(truncnorm)

#Set seed
set.seed(41)

# Helper functions (these need to be defined before parallel processing)
extract_co_max <- function(event_ids, data, window) {
  sapply(event_ids, function(id) {
    start_idx <- max(1, id - window/2)
    end_idx <- min(length(data), id + window/2)
    max(data[start_idx:end_idx], na.rm = TRUE)
  })
}

extract_co_value <- function(event_ids, data, window) {
  sapply(event_ids, function(id) {
    start_idx <- max(1, id - window/2)
    end_idx <- min(length(data), id + window/2)
    max(data[start_idx:end_idx], na.rm = TRUE)
  })
}

extract_co_id <- function(event_ids, data, window) {
  sapply(event_ids, function(id) {
    start_idx <- max(1, id - window/2)
    end_idx <- min(length(data), id + window/2)
    value <- max(data[start_idx:end_idx], na.rm = TRUE)
    matches = which(data[start_idx:end_idx] == value)
    match_select <- if(length(matches) == 1) matches else sample(matches, 1)
    start_idx + match_select - 1
  })
}


classify_event_types <- function(event_ids, names_data, start_idx, end_idx) {
  sapply(seq_along(event_ids), function(i) {
    event_range <- start_idx[i]:end_idx[i]
    ifelse(any(!is.na(names_data[event_range])), "tc", "non_tc")
  })
}


#Function to calculate correlations
calculate_correlations <- function(event_ids, type_ids, precip_data, wl_data) {
  tc_mask <- type_ids == "tc"
  non_tc_mask <- type_ids == "non_tc"
  
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
  non_tc_mask <- type_ids == "non_tc"
  
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
  non_tc_mask <- type_ids == "non_tc"
  
  cop_tc <- if(sum(tc_mask) > 1) {
    tryCatch(
      BiCopSelect(u1 = pobs(precip_data[tc_mask]),
                  u2 = pobs(wl_data[tc_mask])),
      familyset = c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36),
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


# Load data
name = 'Mullica'
year = 1979:2024

#Rainfall event duration
rain_dur = 24

ba = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'/',name,'_aorc_ba_ts.csv',sep=""),header=F)[,1]
ba_lat_lon = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'/',name,'_watershed_centroid.csv',sep=""))

cora_path = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")
cora_path_decl = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_declust_ts.csv',sep="")

#Read in ba aorc data 
all_aorc_time =  seq(as.POSIXct("1979-02-01 00:00:00", tz = "UTC"),as.POSIXct("2022-12-31 23:00:00", tz = "UTC"),by = "hour")

# Remove last two years' worth of data as in original code
aorc_precip <- ba

# Create dataframe and apply HURDAT
aorc_precip_df <- data.frame(all_aorc_time, aorc_precip[1:(length(aorc_precip) - 366*24 - 365*24)])
colnames(aorc_precip_df) <- c("date", "precip")

##Identify individual rainfall events 
# Identify 6-hour dry spells
rain_var<- function(data,sep_crit){
  
  no.rain <- rep(NA, length(data))
  for(k in sep_crit:length(data)){
    no.rain[k] <- ifelse(sum(data[(k-(sep_crit-1)):k]) < 0.1, k, NA)
  }
  no.rain <- na.omit(no.rain)
  
  # Calculate volume between consecutive dry spells (these are rainfall events)
  v_all <- numeric(length(no.rain) - 1)
  for(k in 1:(length(no.rain) - 1)){
    v_all[k] <- sum(data[(no.rain[k] + 1):(no.rain[k+1] - 1)])
  }
  
  # Filter for events with actual rain (v > 0)
  valid_events <- v_all > 0
  
  idx <- which(valid_events)
  idx <- idx[idx + 1 <= length(no.rain)]
  
  start <- no.rain[idx] + 1
  end   <- no.rain[idx + 1] - 1
  dur   <- end - start + 1
  v     <- v_all[idx]
  
  return(list("start" = start,"end" = end, "d" = dur, "v"=v))
}

rainfall_events = rain_var(data=aorc_precip_df$precip,sep_crit=6)
hist(rainfall_events$d, xlab="Event length (hrs)")
length(which(rainfall_events$d>48))/length(rainfall_events$d)
plot(rainfall_events$d,rainfall_events$v)

running_rainfall_24hr <- function(rainfall, window = 24) {
  n <- length(rainfall)
  running_total <- numeric(n) 
  
  for (i in window:n) {
    start <- max(1, i - window +1)
    end   <- i
    running_total[i] <- sum(rainfall[start:end], na.rm = TRUE)
  }
  
  return(running_total)
}

rain_sum_24hr = running_rainfall_24hr(aorc_precip_df$precip)

#rain_sum_24hr 
totals_24hr_window_func <- function(start,end,precip,rain_sum_24hr){
  sum_24    <- numeric(length(start))   
  sum_24_id <- numeric(length(start))  
  hourly_peak <- numeric(length(start))  
  hourly_peak_id <- numeric(length(start))  
  for(i in 1:length(start)){
    sum_24[i] =  max(rain_sum_24hr[start[i]:end[i]])
    sum_24_id[i] = start[i] + which( rain_sum_24hr[start[i]:end[i]] == sum_24[i])[1] - 1
    hourly_peak[i] = max(precip[start[i]:end[i]])
    hourly_peak_id[i] =  start[i] + which(precip[start[i]:end[i]] == hourly_peak[i])[1] - 1
  }
  return(list("sum_24"= sum_24, "sum_24_id"=sum_24_id, "hourly_peak"=hourly_peak, "hourly_peak_id"=hourly_peak_id))
}

totals_24hr_window = totals_24hr_window_func(rainfall_events$start, rainfall_events$end, precip=aorc_precip_df$precip,rain_sum_24hr=rain_sum_24hr)

event_24hr_func <- function(start,end,precip){
  event_sums = numeric(length(precip))
  event_24 = numeric(length(start))
  event_24_id = numeric(length(end))
  for(i in 1:length(start)){
    s = start[i]
    e = end[i]
    len = e-s+1
    
    for(j in 1:len){
      if(j>24){
        event_sums[s+j-1] = sum(precip[(s+j-24):(s+j-1)])
      } else{
        event_sums[s+j-1] = sum(precip[s:(s+j-1)])
      }
    }
    event_24[i] = max(event_sums[s:e])
    event_24_id[i] =  s + which(event_sums[s:e] == event_24[i])[1] - 1
  }
  return(list(event_24=event_24, event_24_id=event_24_id))
}

event_24hr = event_24hr_func(rainfall_events$start, rainfall_events$end, precip=aorc_precip_df$precip)


rain_events = data.frame("start"     = rainfall_events$start, 
                         "end"       = rainfall_events$end, 
                         "dur"       = rainfall_events$d,
                         "vol"       = rainfall_events$v,
                         "sum_24"    = totals_24hr_window$sum_24, 
                         "sum_24_id" = totals_24hr_window$sum_24_id,
                         "hourly_peak"= totals_24hr_window$hourly_peak, 
                         "hourly_peak_id"= totals_24hr_window$hourly_peak_id,
                         "event_24" = event_24hr$event_24,
                         "event_24_id"= event_24hr$event_24_id)

#Decluster rainfall volume
aorc_precip_df <- HURDAT(aorc_precip_df,
                         lat.loc = ba_lat_lon$lat,
                         lon.loc = -ba_lat_lon$lon,
                         rad = 350)

# Decluster events
aorc_precip_df_decl <- Decluster_S_SW(Data = aorc_precip_df[,1:2],
                                      Window_Width_Sum = rain_dur,
                                      Window_Width = 5*24)

#Read in declustered water levels
cora_ntr_decl_df = read.csv(cora_path_decl)[,-1][,c(4,5)]
colnames(cora_ntr_decl_df) <- c("date", "ntr")
cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)

cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                           lat.loc = ba_lat_lon$lat,
                           lon.loc = -ba_lat_lon$lon,
                           rad = 350)

#Find threshold of 24hr totals and ntr peaks that gives 800 events for time lags

find_threshold <- function(decl_data, events_data, target=800) {
  
  for(i in seq(1,0,-0.001)) {
    
    q <- quantile(decl_data, i,  na.rm=TRUE)
    len <- length(which(events_data > q))
    
    if(len >= target) {
      cat("Quantile:", i, "\n")
      cat("Threshold:", q, "\n")
      cat("Number of events:", len, "\n")
      return(q)
    }
  }
  warning("Could not find threshold giving 800 events")
}

rain_threshold_lag <- find_threshold(aorc_precip_df_decl$Totals, rain_events$sum_24)
ntr_threshold_lag <- find_threshold(cora_ntr_decl_df$ntr, cora_ntr_decl_df$ntr)


joint_rp_func <- function(ba, ba_lat_lon, years, cora_path, cora_path_decl, name, rp_rain, rp_ntr, rain_dur, rain_events, rain_thres_lag, ntr_thres_lag, N){
  
  
  #Read in ba aorc data 
  all_aorc_time =  seq(as.POSIXct("1979-02-01 00:00:00", tz = "UTC"),as.POSIXct("2022-12-31 23:00:00", tz = "UTC"),by = "hour")
  
  # Remove last two years' worth of data as in original code
  aorc_precip <- ba
  
  # Create dataframe and apply HURDAT
  aorc_precip_df <- data.frame(all_aorc_time, aorc_precip[1:(length(aorc_precip) - 366*24 - 365*24)])
  colnames(aorc_precip_df) <- c("date", "precip")
  
  
  #Decluster rainfall volume
  aorc_precip_df <- HURDAT(aorc_precip_df,
                           lat.loc = ba_lat_lon$lat,
                           lon.loc = -ba_lat_lon$lon,
                           rad = 350)
  
  # Decluster events
  aorc_precip_df_decl <- Decluster_S_SW(Data = aorc_precip_df[,1:2],
                                        Window_Width_Sum = rain_dur,
                                        Window_Width = 5*24)
  
  # Event identification
  ids <- order(aorc_precip_df_decl$Declustered, decreasing = TRUE)[1:(5*44)]
  types <- classify_event_types(ids, aorc_precip_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(aorc_precip_df), ids + 24))
  rain_event_data = list(ids = ids, types = types)
  
  #Fit GPD
  gpd_aorc = GPD_Fit(Data=aorc_precip_df_decl$Declustered[rain_event_data$ids], Data_Full=aorc_precip_df_decl$Declustered, Thres=min(aorc_precip_df_decl$Declustered[rain_event_data$ids]), Method = "Solari")
  gpd_aorc$Rate = length(rain_event_data$ids)/(nrow(aorc_precip_df)/(24*365.25))
  
  #Read water level time series
  cora_ntr_df = read.csv(cora_path)[,c(1,4)]
  colnames(cora_ntr_df) <- c("date", "ntr")
  cora_ntr_df$date = as.POSIXct(round((cora_ntr_df$date - 719529) * 86400), 
                                origin = "1970-01-01", tz = "UTC")
  
  
  #Read in declustered water levels
  cora_ntr_decl_df = read.csv(cora_path_decl)[,-1][,c(4,5)]
  colnames(cora_ntr_decl_df) <- c("date", "ntr")
  cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)
  
  
  cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                             lat.loc = ba_lat_lon$lat,
                             lon.loc = -ba_lat_lon$lon,
                             rad = 350)
  
  # Event identification
  
  ids <- order(cora_ntr_decl_df$ntr, decreasing = TRUE)[1:(5*44)]
  types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 1*24))
  ntr_event_data = list(ids = ids, types = types)
  
  evt = ntr_event_data
  ntr_ext = evt$ids
  
  #Fit GPD
  gpd_ntr = GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext]), Method="Solari")
  gpd_ntr$Rate = length(cora_ntr_decl_df$ntr[ntr_ext])/(nrow(cora_ntr_df)/(24*365.25))
  
  v_ntr <- extract_co_max(ntr_ext, aorc_precip_df_decl$Totals, window = 3*24)
  v_ntr_id <- extract_co_id(ntr_ext, aorc_precip_df_decl$Totals, window = 3*24)
  
  correlations <- calculate_correlations(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])
  
  cor_ntr_tc <- as.numeric(correlations$tc)
  cor_ntr_non_tc <- as.numeric(correlations$non_tc)
  
  copula = calculate_copula(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])
  
  cop_ntr_tc_family <- as.numeric(copula$tc[1])
  cop_ntr_non_tc_family <- as.numeric(copula$non_tc[1])
  
  cop_ntr_tc_par1 <- as.numeric(copula$tc[2])
  cop_ntr_non_tc_par1 <- as.numeric(copula$non_tc[2])
  
  cop_ntr_tc_par2 <- as.numeric(copula$tc[3])
  cop_ntr_non_tc_par2 <- as.numeric(copula$non_tc[3])
  
  tc_indices <- evt$types == "tc"
  non_tc_indices <- evt$types == "non_tc"
  
  n_ntr_tc = sum(tc_indices)
  n_ntr_non_tc = sum(non_tc_indices)
  
  cop_ntr_tc =BiCopSelect(u1 = pobs(v_ntr[tc_indices]),
                          u2 = pobs(cora_ntr_df$ntr[ntr_ext][tc_indices]),
                          familyset = c(1, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36))
  
  cop_ntr_non_tc =BiCopSelect(u1 = pobs(v_ntr[non_tc_indices]),
                              u2 = pobs(cora_ntr_df$ntr[ntr_ext][non_tc_indices]))
  
  #Fit GPDs
  gpd_ntr_tc = GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext][tc_indices], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext][tc_indices]), Method = "Solari", PLOT = FALSE)
  gpd_ntr_non_tc = GPD_Fit(Data=cora_ntr_decl_df$ntr[ntr_ext][non_tc_indices], Data_Full=cora_ntr_df$ntr, Thres=min(cora_ntr_decl_df$ntr[ntr_ext][non_tc_indices]), Method = "Solari",  PLOT = FALSE)
  
  gpd_ntr_tc$Rate =  n_ntr_tc/44
  gpd_ntr_non_tc$Rate =  n_ntr_non_tc/44
  
  v_ntr_tc = v_ntr[tc_indices] + runif(length(v_ntr[tc_indices]),0.0001,0.001)
  v_ntr_non_tc = v_ntr[non_tc_indices] + runif(length(v_ntr[non_tc_indices]),0.0001,0.001)
  
  #Find best fitting distribution for the conditiionec variable
  non_con_dist_r_tc = Diag_Non_Con_Trunc_AIC(v_ntr_tc, Omit=c("Twe","Weib"))$Best_fit
  non_con_dist_r_non_tc = Diag_Non_Con_Trunc_AIC(v_ntr_non_tc, Omit=c("Twe","Weib"))$Best_fit
  
  #Find non-exceedance probability
  #Obtaining sample conditioned on rainfall
  evt <- rain_event_data
  ntr_vals <- extract_co_max(evt$ids, cora_ntr_df$ntr, window=3*24)
  ntr_vals_id <-extract_co_id(evt$ids, cora_ntr_df$ntr, window=3*24)
  
  samples <-data.frame( "ids_con_rain" = rain_event_data$ids,
                        "type_con_rain" = rain_event_data$types,
                        "rain_con_rain" = aorc_precip_df_decl$Totals[rain_event_data$ids],
                        "ntr_con_rain" = ntr_vals,
                        "ids_ntr_con_rain" = ntr_vals_id,
                        "ids_con_ntr" = ntr_event_data$ids,
                        "type_con_ntr" = ntr_event_data$types,
                        "rain_con_ntr" = v_ntr,
                        "ntr_con_ntr" = cora_ntr_df$ntr[ntr_event_data$ids],
                        "ids_rain_con_ntr" = v_ntr_id)
  
  cat("ntr_vals:", ntr_vals, "\n")
  correlations <- calculate_correlations(evt$ids, evt$types,
                                         aorc_precip_df_decl$Declustered[evt$ids],
                                         ntr_vals)
  
  cor_r_tc <- as.numeric(correlations$tc)
  cor_r_non_tc <- as.numeric(correlations$non_tc)
  
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
  non_tc_r_indices <- evt$types == "non_tc"
  
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
  gpd_rain_tc$Rate = sum(tc_r_indices)/44
  gpd_rain_non_tc$Rate =  sum(non_tc_r_indices)/44
  

  non_con_dist_ntr_tc = Diag_Non_Con_AIC(ntr_vals[tc_r_indices],Omit="Lapl")$Best_fit
  non_con_dist_ntr_non_tc = Diag_Non_Con_AIC(ntr_vals[non_tc_r_indices],Omit="Lapl")$Best_fit
  
  
  #Simulations
  #Number of observations from each copula
  n_sim_r_tc= round(N * (n_r_tc / (n_r_tc + n_ntr_tc + n_r_non_tc + n_ntr_non_tc)),0)
  n_sim_r_non_tc= round(N * (n_r_non_tc / (n_r_tc + n_ntr_tc + n_r_non_tc + n_ntr_non_tc)),0)
  n_sim_ntr_tc= round(N * (n_ntr_tc / (n_r_tc + n_ntr_tc + n_r_non_tc + n_ntr_non_tc)),0)
  n_sim_ntr_non_tc= round(N * (n_ntr_non_tc / (n_r_tc + n_ntr_tc + n_r_non_tc + n_ntr_non_tc)),0)
  
  print(n_sim_r_tc) 
  print(cop_r_tc) 
  
  #Simulate from th fitted copulas
  u_sim_r_tc = data.frame(BiCopSim(n_sim_r_tc, obj=cop_r_tc)) 
  u_sim_r_non_tc = data.frame(BiCopSim(n_sim_r_non_tc, obj=cop_r_non_tc)) 
  u_sim_ntr_tc = data.frame(BiCopSim(n_sim_ntr_tc, obj=cop_ntr_tc)) 
  u_sim_ntr_non_tc = data.frame(BiCopSim(n_sim_ntr_non_tc, obj=cop_ntr_non_tc)) 
  
  #Naming columns
  colnames(u_sim_r_tc) = c("precip","ntr")
  colnames(u_sim_r_non_tc) = c("precip","ntr")
  colnames(u_sim_ntr_tc) = c("precip","ntr")
  colnames(u_sim_ntr_non_tc) = c("precip","ntr") 
  
  #Conditioning variables
  if (gpd_rain_tc$xi != 0) {
    r_sim_r_tc <- ((1 - u_sim_r_tc$precip)^(-gpd_rain_tc$xi) - 1) * gpd_rain_tc$sigma/gpd_rain_tc$xi + gpd_rain_tc$Threshold 
  } else {
    r_sim_r_tc <- gpd_rain_tc$Threshold - gpd_rain_tc$sigma * log(1 - u_sim_r_tc$precip)
  }
  print("r_sim_r_tc")
  print(gpd_rain_tc)
  print(summary(r_sim_r_tc))
  if (gpd_rain_non_tc$xi != 0) {
    r_sim_r_non_tc <- ((1 - u_sim_r_non_tc$precip)^(-gpd_rain_non_tc$xi) - 1) * gpd_rain_non_tc$sigma/gpd_rain_non_tc$xi + gpd_rain_non_tc$Threshold 
  } else {
    r_sim_r_non_tc <- gpd_rain_non_tc$Threshold - gpd_rain_non_tc$sigma * log(1 - u_sim_r_non_tc$precip)
  }
  
  print("r_sim_r_non_tc")
  print(summary(r_sim_r_non_tc))
  if (gpd_ntr_tc$xi != 0) {
    ntr_sim_ntr_tc <- ((1 - u_sim_ntr_tc$ntr)^(-gpd_ntr_tc$xi) - 1) * gpd_ntr_tc$sigma/gpd_ntr_tc$xi + gpd_ntr_tc$Threshold 
  } else {
    ntr_sim_ntr_tc <- gpd_ntr_tc$Threshold - gpd_ntr_tc$sigma * log(1 - u_sim_ntr_tc$ntr)
  }
  
  print("ntr_sim_ntr_tc")
  print(summary(ntr_sim_ntr_tc))
  if (gpd_ntr_non_tc$xi != 0) {
    ntr_sim_ntr_non_tc <- ((1 - u_sim_ntr_non_tc$ntr)^(-gpd_ntr_non_tc$xi) - 1) * gpd_ntr_non_tc$sigma/gpd_ntr_non_tc$xi + gpd_ntr_non_tc$Threshold 
  } else {
    ntr_sim_ntr_non_tc <- gpd_ntr_non_tc$Threshold - gpd_ntr_non_tc$sigma * log(1 - u_sim_ntr_non_tc$ntr)
  }
  
  print("ntr_sim_ntr_non_tc")
  print(summary(ntr_sim_ntr_non_tc))
  #Conditioned variables
  #Marginal distributions for ntr in con r samples
  #Tc
  if(non_con_dist_ntr_tc == "Gum"){
    non_con_dist_fit_ntr_tc <- gamlss(ntr_vals[tc_r_indices]  ~ 1, family= GU)
    ntr_sim_r_tc<-qGU(u_sim_r_tc$ntr,as.numeric(non_con_dist_fit_ntr_tc$mu.coefficients),exp(as.numeric(non_con_dist_fit_ntr_tc$sigma.coefficients)))
  }
  if(non_con_dist_ntr_tc=="RGum"){
    non_con_dist_fit_ntr_tc <- gamlss(ntr_vals[tc_r_indices] ~ 1,family=RG)
    ntr_sim_r_tc<-qRG(u_sim_r_tc$ntr,non_con_dist_fit_ntr_tc$mu.coefficients,exp(non_con_dist_fit_ntr_tc$sigma.coefficients))
  }
  if(non_con_dist_ntr_tc=="Gaus"){
    non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices],"normal")
    ntr_sim_r_tc<-qnorm(u_sim_r_tc$ntr, as.numeric(non_con_dist_fit_ntr_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
  }
  if(non_con_dist_ntr_tc=="Lapl"){
    non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices], dlaplace, start=list(location=mean(ntr_vals[tc_r_indices]), scale=sd(ntr_vals[tc_r_indices])/sqrt(2)))
    ntr_sim_r_tc <- qlaplace(u_sim_r_tc$ntr,as.numeric(non_con_dist_fit_ntr_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
  }
  if(non_con_dist_ntr_tc=="Logis"){
    non_con_dist_fit_ntr_tc<-fitdistr(ntr_vals[tc_r_indices],"logistic")
    ntr_sim_r_tc<-qlogis(u_sim_r_tc$ntr,as.numeric(non_con_dist_fit_ntr_tc$estimate[1]),as.numeric(non_con_dist_fit_ntr_tc$estimate[2]))
  }
  
  print("ntr_sim_r_tc")
  print(summary(ntr_sim_r_tc))
  
  #Non-Tc
  if(non_con_dist_ntr_non_tc == "Gum"){
    non_con_dist_fit_ntr_non_tc <- gamlss(ntr_vals[non_tc_r_indices]  ~ 1, family= GU)
    ntr_sim_r_non_tc<-qGU(u_sim_r_non_tc$ntr,as.numeric(non_con_dist_fit_ntr_non_tc$mu.coefficients),exp(as.numeric(non_con_dist_fit_ntr_non_tc$sigma.coefficients)))
  }
  if(non_con_dist_ntr_non_tc=="RGum"){
    non_con_dist_fit_ntr_non_tc <- gamlss(ntr_vals[non_tc_r_indices] ~ 1,family=RG)
    ntr_sim_r_non_tc<-qRG(u_sim_r_non_tc$ntr,non_con_dist_fit_ntr_non_tc$mu.coefficients,exp(non_con_dist_fit_ntr_non_tc$sigma.coefficients))
  }
  if(non_con_dist_ntr_non_tc=="Gaus"){
    non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices],"normal")
    ntr_sim_r_non_tc<-qnorm(u_sim_r_non_tc$ntr, as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
  }
  if(non_con_dist_ntr_non_tc=="Lapl"){
    non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices], dlaplace, start=list(location=mean(ntr_vals[non_tc_r_indices]), scale=sd(ntr_vals[non_tc_r_indices])/sqrt(2)))
    ntr_sim_r_non_tc <- qlaplace(u_sim_r_non_tc$ntr,as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]), as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
  }
  if(non_con_dist_ntr_non_tc=="Logis"){
    non_con_dist_fit_ntr_non_tc<-fitdistr(ntr_vals[non_tc_r_indices],"logistic")
    ntr_sim_r_non_tc<-qlogis(u_sim_r_non_tc$ntr,as.numeric(non_con_dist_fit_ntr_non_tc$estimate[1]),as.numeric(non_con_dist_fit_ntr_non_tc$estimate[2]))
  }
  print("ntr_sim_r_non_tc")
  print(summary(ntr_sim_r_non_tc))
  
  #Marginal distributions for rainfall in sample con. on ntr
  #Tc
  if(non_con_dist_r_tc=="BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = v_ntr_tc)
    non_con_dist_fit_r_tc<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    r_sim_ntr_tc<-qbisa(u_sim_ntr_tc$precip, as.numeric(Coef(non_con_dist_fit_r_tc)[1]),as.numeric(Coef(non_con_dist_fit_r_tc)[2]))
  }
  if(non_con_dist_r_tc=="Exp"){
    non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc,"exponential")
    r_sim_ntr_tc<-qexp(u_sim_ntr_tc$precip, as.numeric(non_con_dist_fit_r_tc$estimate[1]))
  }
  if(non_con_dist_r_tc=="Gam(2)"){
    non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc, "gamma")
    r_sim_ntr_tc<-qgamma(u_sim_ntr_tc$precip, shape = as.numeric(non_con_dist_fit_r_tc$estimate[1]), rate = as.numeric(non_con_dist_fit_r_tc$estimate[2]))
  }
  if(non_con_dist_r_tc=="Gam(3)"){
    data.gamlss = data.frame(X=v_ntr_tc)
    non_con_dist_fit_r_tc <-  tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                       error = function(e) "error")
    r_sim_ntr_tc<-qGG(u_sim_ntr_tc$precip, mu=exp(non_con_dist_fit_r_tc$mu.coefficients), sigma=exp(non_con_dist_fit_r_tc$sigma.coefficients), nu=non_con_dist_fit_r_tc$nu.coefficients)
  }
  if(non_con_dist_r_tc=="InvG"){
    non_con_dist_fit_r_tc <- fitdist(v_ntr_tc, "invgauss", start = list(mean = 5, shape = 1))
    r_sim_ntr_tc<-qinvgauss(u_sim_ntr_tc$precip, as.numeric(non_con_dist_fit_r_tc$estimate[1]), as.numeric(non_con_dist_fit_r_tc$estimate[2]))
  }
  if(non_con_dist_r_tc=="LogN"){
    non_con_dist_fit_r_tc <- fitdistr(v_ntr_tc,"lognormal")
    r_sim_ntr_tc<-qlnorm(u_sim_ntr_tc$precip, meanlog = as.numeric(non_con_dist_fit_r_tc$estimate[1]), sdlog = as.numeric(non_con_dist_fit_r_tc$estimate[2]))
  }
  if(non_con_dist_r_tc=="TNorm"){
    non_con_dist_fit_r_tc <-fitdistr(v_ntr_tc,"normal")
    r_sim_ntr_tc<-qtruncnorm(u_sim_ntr_tc$precip,a=min(v_ntr_tc),as.numeric(non_con_dist_fit_r_tc$estimate[1]),as.numeric(non_con_dist_fit_r_tc$estimate[2]))
  }
  if(non_con_dist_r_tc=="Twe"){
    non_con_dist_fit_r_tc <-tweedie.profile(v_ntr_tc ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    r_sim_ntr_tc<-qtweedie(u_sim_ntr_tc$precip, power=non_con_dist_fit_r_tc$p.max, mu=mean(v_ntr_tc), phi=non_con_dist_fit_r_tc$phi.max)
  }
  if(non_con_dist_r_tc=="Weib"){
    non_con_dist_fit_r_tc <- fitdistr(v_ntr_tc, "weibull")
    r_sim_ntr_tc<-qweibull(u_sim_ntr_tc$precip, as.numeric(non_con_dist_fit_r_tc$estimate[1]), as.numeric(non_con_dist_fit_r_tc$estimate[2]))
  }
  print("r_sim_ntr_tc")
  print(summary(r_sim_ntr_tc))
  
  #non-Tc
  if(non_con_dist_r_non_tc=="BS"){
    bdata2 <- data.frame(shape = exp(-0.5), scale = exp(0.5))
    bdata2 <- transform(bdata2, y = v_ntr_non_tc)
    non_con_dist_fit_r_non_tc<-vglm(y ~ 1, bisa, data = bdata2, trace = FALSE)
    r_sim_ntr_non_tc<-qbisa(u_sim_ntr_non_tc$precip, as.numeric(Coef(non_con_dist_fit_r_non_tc)[1]),as.numeric(Coef(non_con_dist_fit_r_non_tc)[2]))
  }
  if(non_con_dist_r_non_tc=="Exp"){
    non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc,"exponential")
    r_sim_ntr_non_tc<-qexp(u_sim_ntr_non_tc$precip, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]))
  }
  if(non_con_dist_r_non_tc=="Gam(2)"){
    non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc, "gamma")
    r_sim_ntr_non_tc<-qgamma(u_sim_ntr_non_tc$precip, shape = as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), rate = as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
  }
  if(non_con_dist_r_non_tc=="Gam(3)"){
    data.gamlss = data.frame(X=v_ntr_non_tc)
    non_con_dist_fit_r_non_tc <-  tryCatch(gamlss(X~1, data=data.gamlss, family=GG),
                                           error = function(e) "error")
    r_sim_ntr_non_tc<-qGG(u_sim_ntr_non_tc$precip, mu=exp(non_con_dist_fit_r_non_tc$mu.coefficients), sigma=exp(non_con_dist_fit_r_non_tc$sigma.coefficients), nu=non_con_dist_fit_r_non_tc$nu.coefficients)
  }
  if(non_con_dist_r_non_tc=="InvG"){
    non_con_dist_fit_r_non_tc <- fitdist(v_ntr_non_tc, "invgauss", start = list(mean = 5, shape = 1))
    r_sim_ntr_non_tc<-qinvgauss(u_sim_ntr_non_tc$precip, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
  }
  if(non_con_dist_r_non_tc=="LogN"){
    non_con_dist_fit_r_non_tc <- fitdistr(v_ntr_non_tc,"lognormal")
    r_sim_ntr_non_tc<-qlnorm(u_sim_ntr_non_tc$precip, meanlog = as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), sdlog = as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
  }
  if(non_con_dist_r_non_tc=="TNorm"){
    non_con_dist_fit_r_non_tc <-fitdistr(v_ntr_non_tc,"normal")
    r_sim_ntr_non_tc<-qtruncnorm(u_sim_ntr_non_tc$precip,a=min(v_ntr_non_tc),as.numeric(non_con_dist_fit_r_non_tc$estimate[1]),as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
  }
  if(non_con_dist_r_non_tc=="Twe"){
    non_con_dist_fit_r_non_tc <-tweedie.profile(v_ntr_non_tc ~ 1,p.vec=seq(1.5, 2.5, by=0.2), do.plot=FALSE)
    r_sim_ntr_non_tc<-qtweedie(u_sim_ntr_non_tc$precip, power=non_con_dist_fit_r_non_tc$p.max, mu=mean(v_ntr_non_tc), phi=non_con_dist_fit_r_non_tc$phi.max)
  }
  if(non_con_dist_r_non_tc=="Weib"){
    non_con_dist_fit_r_non_tc <- fitdistr(v_ntr_non_tc, "weibull")
    r_sim_ntr_non_tc<-qweibull(u_sim_ntr_non_tc$precip, as.numeric(non_con_dist_fit_r_non_tc$estimate[1]), as.numeric(non_con_dist_fit_r_non_tc$estimate[2]))
  }

  #Extracting rainfall and ntr lag times 
  rain_above = which(rain_events$event_24>rain_thres_lag)
  rain_events_filt = rain_events[rain_above, ]

  rain_minus_ntr_rain  = extract_co_id(event_ids=rain_events_filt$hourly_peak_id, data = cora_ntr_df$ntr, window=3*24)
  rain_minus_ntr_rain_ntr_val  = extract_co_value(event_ids=rain_events_filt$hourly_peak_id, data = cora_ntr_df$ntr, window=3*24)

  rain_minus_ntr_rain = rain_events_filt$hourly_peak_id - rain_minus_ntr_rain 

  
  #Only keep ntrs above 50% percentile
  cora_ntr_decl_df_EventID = which(!is.na(cora_ntr_decl_df$ntr))
  ntr_above = which(cora_ntr_decl_df$ntr[cora_ntr_decl_df_EventID]>ntr_thres_lag)
  cora_ntr_decl_df_EventID = cora_ntr_decl_df_EventID[ntr_above]
  cora_ntr_decl_df_EventID = cora_ntr_decl_df_EventID[cora_ntr_decl_df_EventID>110]
  cora_ntr_decl_df_EventID = cora_ntr_decl_df_EventID[cora_ntr_decl_df_EventID<(length(cora_ntr_decl_df_EventID)-100)]
  rain_minus_ntr_ntr = extract_co_id(event_ids=cora_ntr_decl_df_EventID, data = aorc_precip_df$precip, window=3*24)
  rain_minus_ntr_ntr = rain_minus_ntr_ntr - cora_ntr_decl_df_EventID

  
  #print("rain_minus_ntr_ntr")
  #print(summary(rain_minus_ntr_ntr))
  
  #Putting simulations into a dataframe
  sim = data.frame(c(r_sim_r_tc,r_sim_r_non_tc,r_sim_ntr_tc,r_sim_ntr_non_tc),c(ntr_sim_r_tc,ntr_sim_r_non_tc,ntr_sim_ntr_tc,ntr_sim_ntr_non_tc),c(rep("r_tc",n_sim_r_tc),rep("r_non_tc",n_sim_r_non_tc),rep("ntr_tc",n_sim_ntr_tc),rep("ntr_non_tc",n_sim_ntr_non_tc)))
  colnames(sim) = c("precip","ntr","sample")
  
  #Generating rainfall fields
  rain_event_id = numeric(nrow(sim))
  rain_factor = numeric(nrow(sim))
  rain_minus_ntr_rain_event = numeric(nrow(sim))
  print("Here_1")
  for(m in 1:nrow(sim)){
    # Initialize rain_factor for this iteration
    rain_factor[m] <- Inf  # Start with a value >= 4
    #while(rain_factor[m]>4){
    #probabilities of sampling an observed event as inverse of distance between observed 24-hr rainfall of delclustered events and simulated 24-hr rainfall
    d <- abs(sim$precip[m]-rain_events_filt$event_24)
    w <- 1 / (d + 1e-10)
    probs <- w / sum(w)
    row_idx <- sample(1:nrow(rain_events_filt), size=1, replace=TRUE, prob=probs)
    rain_event_id[m] <- rain_events_filt$event_24_id[row_idx]  # store the timestep
    rain_factor[m] <- sim$precip[m] / rain_events_filt$event_24[row_idx] 
    rain_minus_ntr_rain_event[m] = rain_minus_ntr_rain[row_idx]
    #}
  }
  
  #Generating NTR peak
  ntr_event_id = numeric(nrow(sim))
  ntr_factor = numeric(nrow(sim))
  rain_minus_ntr_ntr_event = numeric(nrow(sim))
  print("here_2")
  n_reject = numeric(nrow(sim))
  for(n in 1:nrow(sim)){
    
    #probabilities of sampling an observed event as inverse of distance between observed (declustered) ntr and simulated peak ntr
    d <-  abs(sim$ntr[n] - cora_ntr_decl_df$ntr[cora_ntr_decl_df_EventID])
    w <- 1 / (d + 1e-10)
    probs <- w / sum(w)
    
    no_secondary_peak = FALSE
    n_reject[n] = 0
    while(!no_secondary_peak){
      
      # Sample a position within cora_ntr_decl_df_EventID
      row_idx <- sample(length(cora_ntr_decl_df_EventID), size=1, replace=TRUE, prob=probs)
      ntr_event_id[n] <- cora_ntr_decl_df_EventID[row_idx]  # store the timestep
      ntr_factor[n] <-sim$ntr[n] / cora_ntr_decl_df$ntr[ntr_event_id[n]]
      rain_minus_ntr_ntr_event[n] <- rain_minus_ntr_ntr[row_idx]  # index by position
      
      if(sim$sample[n]=="r_tc" | sim$sample[n]=="r_non_tc"){
        lead <- max(1, ntr_event_id[n] - 72 + rain_minus_ntr_ntr_event[n])
        tail <- min(length(cora_ntr_df), ntr_event_id[n] + 72 + rain_minus_ntr_ntr_event[n])
      } else{
        lead <- max(1,                   ntr_event_id[n] - 72)
        tail <- min(length(cora_ntr_df), ntr_event_id[n] + 72)
      }
      
      no_secondary_peak <- max(cora_ntr_df$ntr[lead:tail]) <= cora_ntr_df$ntr[ntr_event_id[n]] + 1e-10
      if(!no_secondary_peak) n_reject[n] = n_reject[n] + 1
    }  
  }
  #Putting simulations into a dataframe
  sim = cbind(sim,"rain_id" = rain_event_id, "rain_factor" = rain_factor, "rain_minus_ntr_rain_event"=rain_minus_ntr_rain_event, 
              "ntr_id" = ntr_event_id, "ntr_factor" = ntr_factor, "rain_minus_ntr_ntr_event"=rain_minus_ntr_ntr_event, "n_reject"=n_reject)
  
  
  # Store results
  res <- list("res" = data.frame(
    cora = 3,
    cor_r_tc = cor_r_tc, cor_r_non_tc = cor_r_non_tc,
    cor_ntr_tc = cor_ntr_tc, cor_ntr_non_tc = cor_ntr_non_tc,
    cop_r_tc_family = cop_r_tc_family, cop_r_non_tc_family = cop_r_non_tc_family,
    cop_ntr_tc_family  = cop_ntr_tc_family, cop_ntr_non_tc_family = cop_ntr_non_tc_family,
    cop_r_tc_par1 = cop_r_tc_par1 , cop_r_non_tc_par1 = cop_r_non_tc_par1,
    cop_ntr_tc_par1  = cop_ntr_tc_par1,cop_ntr_non_tc_par1 = cop_ntr_non_tc_par1,
    cop_r_tc_par2 = cop_r_tc_par2, cop_r_non_tc_par2 = cop_r_non_tc_par2,
    cop_ntr_tc_par2  = cop_ntr_tc_par2, cop_ntr_non_tc_par2 = cop_ntr_non_tc_par2,
    n_r_tc = n_r_tc, n_r_non_tc = n_r_non_tc,
    n_ntr_tc = n_ntr_tc, n_ntr_non_tc  = n_ntr_non_tc,
    non_con_dist_r_tc = non_con_dist_r_tc, non_con_dist_r_non_tc = non_con_dist_r_non_tc,
    non_con_dist_ntr_tc  = non_con_dist_ntr_tc, non_con_dist_ntr_non_tc = non_con_dist_ntr_non_tc),
    "samples" = samples, 
    "rain_minus_ntr_rain" = rain_minus_ntr_rain, 
    "rain_minus_ntr_ntr" = rain_minus_ntr_ntr, 
    "rain_minus_ntr_rain_ntr_val" =rain_minus_ntr_rain_ntr_val,
    "sim" = sim)
  
  return(res)
  
}

res = joint_rp_func(ba=ba, ba_lat_lon= ba_lat_lon, years=years, cora_path=cora_path, cora_path_decl=cora_path_decl, name=name, rp_rain=10, rp_ntr=10, rain_dur=24, rain_events=rain_events, rain_thres_lag=rain_threshold_lag, ntr_thres_lag=ntr_threshold_lag, N=10000)


write.csv(res$res,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res.csv',sep=""))
write.csv(res$samples,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_samples.csv',sep=""))
write.csv(res$rain_minus_ntr_rain,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain.csv',sep=""))
write.csv(res$rain_minus_ntr_ntr,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_ntr.csv',sep=""))
write.csv(res$sim,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim.csv',sep=""))
write.csv(res$rain_minus_ntr_rain_ntr_val,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain_ntr_val.csv',sep=""))

detach("package:MultiHazard", unload = TRUE)
detach("package:VGAM", unload = TRUE)

library('MultiHazard')
library(texmex)

#write.csv(res$res,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res.csv',sep=""))
df_sim = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim.csv',sep=""))
df_samples = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_samples.csv',sep=""))
rain_minus_ntr_rain = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain.csv',sep=""))[,-1]
rain_minus_ntr_ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_ntr.csv',sep=""))[,-1]
rain_minus_ntr_rain_ntr_val = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain_ntr_val.csv',sep=""))[,-1]
summary(rain_minus_ntr_rain_ntr_val)


#rain_minus_ntr = c(df_sample$ids_con_rain - df_sample$ids_ntr_con_rain, df_sample$ids_rain_con_ntr - df_sample$ids_con_ntr)

cora.file = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")
cora = read.csv(cora.file)

cora$date_time <- as.POSIXct(round((cora$date_time - 719529) * 86400 / 3600) * 3600, 
                             origin = "1970-01-01", 
                             tz = "UTC")

#Identify ntr month
#Read water level time series
cora_ntr_df = read.csv(cora_path)[,c(1,4)]
colnames(cora_ntr_df) <- c("date", "ntr")
cora_ntr_df$date = as.POSIXct(cora_ntr_df$date)

#Read in declustered water levels
cora_ntr_decl_df = read.csv(cora_path_decl)[,-1][,c(4,5)]
colnames(cora_ntr_decl_df) <- c("date", "ntr")
cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)

cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                           lat.loc = ba_lat_lon$lat,
                           lon.loc = -ba_lat_lon$lon,
                           rad = 350)

# Event identification

ids <- order(cora_ntr_decl_df$ntr, decreasing = TRUE)[1:(5*44)]
types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 1*24))
ntr_event_data = list(ids = ids, types = types)

cora_df = read.csv(cora_path)#[,c(1,4)]
#colnames(cora_df) <- c("date", "ntr")
cora_df$date_time <- as.POSIXct(round((cora_df$date_time - 719529) * 86400 / 3600) * 3600, 
                                origin = "1970-01-01", 
                                tz = "UTC")

plot((ids[1]-24):(ids[1]+24),cora_df$ntr[(ids[1]-24):(ids[1]+24)], col=1)
points((ids[1]-24):(ids[1]+24),cora_df$pred[(ids[1]-24):(ids[1]+24)], col=2)

##Sampling
#Simulate rain - ntr peak lag

rain_minus_ntr_rain = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain.csv',sep=""))[,-1]
rain_minus_ntr_ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_ntr.csv',sep=""))[,-1]

sim_rain_minus_ntr =rep(NA, nrow(df_sim))
for(i in 1:nrow(df_sim)){
  if(df_sim$sample[i] %in% c("r_tc", "r_non_tc")){
    sim_rain_minus_ntr[i] = df_sim$rain_minus_ntr_rain_event[i] #sample(rain_minus_ntr_rain,size=1)
  } else{
    sim_rain_minus_ntr[i] =  df_sim$rain_minus_ntr_ntr_event[i] #sample(rain_minus_ntr_ntr,size=1)
  }
}

#Simulating month and ntr - tide lag
ntr_count_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_count_tc.csv',sep=""))
ntr_count_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_count_non_tc.csv',sep=""))
ntr_minus_tide_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_minus_tide__tc.csv',sep=""))
ntr_minus_tide_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_minus_tide_non_tc.csv',sep=""))

#Simulate rain - ntr peak lag

rain_minus_ntr_rain = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_rain.csv',sep=""))[,-1]
rain_minus_ntr_ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_minus_ntr_ntr.csv',sep=""))[,-1]

sim_rain_minus_ntr =rep(NA, nrow(df_sim))
for(i in 1:nrow(df_sim)){
  if(df_sim$sample[i] %in% c("r_tc", "r_non_tc")){
    sim_rain_minus_ntr[i] = df_sim$rain_minus_ntr_rain_event[i] #sample(rain_minus_ntr_rain,size=1)
  } else{
    sim_rain_minus_ntr[i] =  df_sim$rain_minus_ntr_ntr_event[i] #sample(rain_minus_ntr_ntr,size=1)
  }
}

#Simulating month and ntr - tide lag
ntr_count_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_count_tc.csv',sep=""))
ntr_count_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_count_non_tc.csv',sep=""))
ntr_minus_tide_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_minus_tide__tc.csv',sep=""))
ntr_minus_tide_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_minus_tide_non_tc.csv',sep=""))
msl_month = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_msl_month_year.csv',sep=""))[-1,]
tide_month = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_tide_month.csv',sep=""))

probs_tc = ntr_count_tc$count/sum(ntr_count_tc$count)
probs_non_tc = ntr_count_non_tc$count/sum(ntr_count_non_tc$count)

sim_month = rep(NA, nrow(df_sim))
sim_msl =  rep(NA, nrow(df_sim))
sim_ntr_minus_tide =rep(NA, nrow(df_sim))
for(i in 1:nrow(df_sim)){
  if(df_sim$sample[i] %in% c("r_tc", "ntr_tc")){
    sim_month[i] = sample(ntr_count_tc$month,size=1,prob=probs_tc)
    sim_ntr_minus_tide[i] = sample(ntr_minus_tide_tc[,1],size=1)
  } else{
    sim_month[i] = sample(ntr_count_non_tc$month,size=1,prob=probs_non_tc)
    sim_ntr_minus_tide[i] =  sample(ntr_minus_tide_non_tc[,1],size=1)
  }
}

plot(df_samples$ntr_con_ntr[c(which(df_samples$type_con_ntr=="tc"),which(df_samples$type_con_ntr=="non_tc"))],
     c(ntr_minus_tide_tc[,1],ntr_minus_tide_non_tc[,1]),
     col=c(rep("Red",length(which(df_samples$type_con_ntr=="tc"))),rep("Black",length(which(df_samples$type_con_ntr=="non_tc")))),
     pch=16,
     xlab="ntr", ylab="Lag (ntr - tide) (hr)")


plot(df_samples$ntr_con_ntr[c(which(df_samples$type_con_ntr=="tc"),which(df_samples$type_con_ntr=="non_tc"))],
     c(ntr_minus_tide_tc[,1],ntr_minus_tide_non_tc[,1]),
     col=c(rep("Red",length(which(df_samples$type_con_ntr=="tc"))),rep("Black",length(which(df_samples$type_con_ntr=="non_tc")))),
     pch=16,
     xlab="ntr", ylab="Lag (ntr - tide) (hr)")

#write.csv(res$res,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_minus_tide.csv',sep=""))
#write.csv(res$samples,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_month.csv',sep=""))

#Which high tides are not in the past 18.6 years
nrow(cora_df) - 18.6*365.25*24

#Remove high tides with elements less than 221912
tide_month = tide_month[-(1:(max(which((tide_month$high_tides<221912))))),]
tide_month = tide_month[-which(tide_month$high_tides>384919),]
tide_month = tide_month[-nrow(tide_month),]
  
sim_tide = rep(NA, nrow(df_sim))
sim_msl = rep(NA, nrow(df_sim))

for(i in 1:nrow(df_sim)){
  sim_msl[i] = sample(msl_month$msl_month[which(msl_month$month == sim_month[i])],size=1)
  sim_tide[i] = sample(tide_month$high_tides[which(tide_month$month == sim_month[i])],size=1)
}


#Generating water level time series
wl_generation <- function(df_sim, sim_tide, sim_ntr_minus_tide,
                          sim_msl, cora_df) {
  N <- nrow(df_sim)
  
  # ── Index vectors ─────────────────────────────────────────────────────────
  pred_center <- as.integer(sim_tide + sim_ntr_minus_tide)   # centre in cora_df$pred
  ntr_center  <- as.integer(df_sim$ntr_id)                   # centre in cora_df$ntr
  
  # ── Full 145-column extraction (vectorised, no loop) ──────────────────────
  # outer() builds an N x 145 index matrix; single vector lookup 
  pred_cols <- outer(pred_center, -72L:72L, `+`)   # N x 145
  ntr_cols  <- outer(ntr_center,  -72L:72L, `+`)   # N x 145
  
  sim_pred <- matrix(cora_df$pred[pred_cols], nrow = N)
  sim_ntr  <- sweep(                                # sweep avoids column-major
    matrix(cora_df$ntr[ntr_cols], nrow = N),        # recycling bug (ntr_factor
    1, df_sim$ntr_factor, `*`)                       # is per-row, not per-col)
  
  # sim_msl is also per-row — same sweep pattern
  sim_wl <- sweep(sim_pred + sim_ntr, 1, sim_msl, `+`)
  
  # ── Trim to 73 columns centred on the appropriate peak ────────────────────
  # rain-conditioned rows shift centre by rain_minus_ntr_rain_event; others stay at 73
  is_rain     <- df_sim$sample %in% c("r_tc", "r_non_tc")
  offset      <- ifelse(is_rain,
                        as.integer(df_sim$rain_minus_ntr_rain_event), 0L)
  trim_center <- 73L + offset                        # column index within 145-col arrays
  
  trim_cols <- outer(trim_center, -36L:36L, `+`)    # N x 73
  
  # Row-wise slice via a two-column integer index matrix
  rc_idx <- cbind(
    rep(seq_len(N), 73L),    # row indices (repeated for each of 73 columns)
    as.integer(trim_cols)    # column indices
  )
  
  extract_trim <- function(mat) matrix(mat[rc_idx], nrow = N, ncol = 73L)
  
  list(
    sim_pred_trim = extract_trim(sim_pred),
    sim_ntr_trim  = extract_trim(sim_ntr),
    sim_wl_trim   = extract_trim(sim_wl)
  )
}

wl = wl_generation(df_sim = df_sim, sim_tide = sim_tide, sim_ntr_minus_tide = sim_ntr_minus_tide,
              sim_msl = sim_msl, cora_df= cora_df)


rainfall_field = function(df_sim, rain_events, aorc_precip,sim_rain_minus_ntr = NULL){
  sim_rain = array(0, dim=c(nrow(df_sim),145))
  for(i in 1:nrow(df_sim)){
    idx = which(rain_events$event_24_id==df_sim$rain_id[i])
    if(df_sim$sample[i]=="r_tc" | df_sim$sample[i]=="r_non_tc"){
      
      lead = max(0,72-(rain_events$hourly_peak_id[idx]-rain_events$start[idx]))
      trail = max(0,72-(rain_events$end[idx]-rain_events$hourly_peak_id[idx]))
      
      rain = df_sim$rain_factor[i] * aorc_precip_df$precip[(rain_events$start[idx]):(rain_events$end[idx])]
      
      if( (rain_events$hourly_peak_id[idx]-rain_events$start[idx])>72){
        n_remove = rain_events$hourly_peak_id[idx]-rain_events$start[idx]-72
        rain = rain[-(1:n_remove)]
      }
      
      if( (rain_events$end[idx]-rain_events$hourly_peak_id[idx])>72){
        n_remove = rain_events$end[idx]-rain_events$hourly_peak_id[idx]-72
        rain = rain[-((length(rain)-n_remove+1):length(rain))]
      }
      
      sim_rain[i,] = c(rep(0,lead),rain,rep(0,trail))
    } else {
      
      
      lead = max(0,(72+sim_rain_minus_ntr[i]-(rain_events$hourly_peak_id[idx]-rain_events$start[idx])))
      tail = max(0,(72-sim_rain_minus_ntr[i]-(rain_events$end[idx]-rain_events$hourly_peak_id[idx])))
      rain = df_sim$rain_factor[i] *aorc_precip_df$precip[(rain_events$start[idx]):(rain_events$end[idx])]
      
      
      if(72+sim_rain_minus_ntr[i] < (rain_events$hourly_peak_id[idx]-rain_events$start[idx])){
        n_remove = (rain_events$hourly_peak_id[idx]-rain_events$start[idx]) - (72+sim_rain_minus_ntr[i]) 
        rain = rain[-(1:n_remove)]
      }
      
      if(72-sim_rain_minus_ntr[i] < (rain_events$end[idx]-rain_events$hourly_peak_id[idx])){
        n_remove = (rain_events$end[idx]-rain_events$hourly_peak_id[idx]) - (72-sim_rain_minus_ntr[i])
        rain = rain[-((length(rain)-n_remove+1):length(rain))]
      }
      
      sim_rain[i,] = c(rep(0,lead),rain,rep(0,tail))
    }
  }
  return(sim_rain)
}



##Rainfall field at AORC subset points

library('ncdf4')

res_rain_list = list()
res_vol_list = list()
res_cor_list = list()

res_rain_df_tc_list = list()
res_vol_df_tc_list = list()

res_rain_df_non_tc_list = list()
res_vol_df_non_tc_list = list()

res_cor_con_rain_list = list()
res_cor_con_ntr_list = list()

res_ddf_list_tc <- list()
res_ddf_list_non_tc <- list()
res_ddf_all_con_r_list <-list()

res_cor_rain_wl_con_rain_list = list()
res_cor_rain_wl_con_wl_list = list()
res_cor_24hr_rain_wl_con_24hr_rain_list = list()
res_cor_24hr_rain_wl_con_wl_list = list()
res_cor_24hr_rain_ntr_con_24hr_rain_list = list()
res_cor_24hr_rain_ntr_con_ntr_list = list()


start <- proc.time()
for(i in 1:2355){
print(i)

# AORC data processing
aorc_precip <- numeric()
aorc_time <- list()


years = 1979:2024

for(j in seq_along(years)) {
  
  aorc_file <- nc_open(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/Mullica/aorc_subset/',name,'_AORC_sub_', years[j], '.nc',sep=""))
  
  aorc_precip_get <- ncvar_get(aorc_file, "APCP_surface_1", start = c(i, 1), 
                               count = c(1, -1))
  
  aorc_precip_get <- aorc_precip_get#[-1]
  aorc_precip <- c(aorc_precip, aorc_precip_get)
  
  time_raw <- ncvar_get(aorc_file, "time")
  
  # Step 3: Convert to POSIXct
  aorc_time_get <- as.POSIXct(time_raw, origin = "1970-01-01", tz = "UTC")
  
  aorc_time[[j]] <- aorc_time_get
  nc_close(aorc_file) 
}
all_aorc_time = do.call(c,aorc_time)

# Remove last two years' worth of data as in original code
all_aorc_time <- all_aorc_time[1:(length(all_aorc_time)-366*24-365*24)]

# Create dataframe and apply HURDAT
aorc_precip_df <- data.frame(all_aorc_time, aorc_precip[1:(length(aorc_precip) - 366*24 - 365*24)])
colnames(aorc_precip_df) <- c("date", "precip")


#Read water level time series
cora_ntr_df = read.csv(cora_path)[,c(1,4)]
colnames(cora_ntr_df) <- c("date", "ntr")
cora_ntr_df$date = as.POSIXct(round((cora_ntr_df$date - 719529) * 86400), 
                              origin = "1970-01-01", tz = "UTC")


##Identify individual rainfall events at AORC point 1
# Identify 6-hour dry spells
rain_var<- function(data,sep_crit){
  
  no.rain <- rep(NA, length(data))
  for(k in sep_crit:length(data)){
    no.rain[k] <- ifelse(sum(data[(k-(sep_crit-1)):k]) == 0, k, NA)
  }
  no.rain <- na.omit(no.rain)
  
  # Calculate volume between consecutive dry spells (these are rainfall events)
  v_all <- numeric(length(no.rain) - 1)
  for(k in 1:(length(no.rain) - 1)){
    v_all[k] <- sum(data[(no.rain[k] + 1):(no.rain[k+1] - 1)])
  }
  
  # Filter for events with actual rain (v > 0)
  valid_events <- v_all > 0
  
  idx <- which(valid_events)
  idx <- idx[idx + 1 <= length(no.rain)]
  
  start <- no.rain[idx] + 1
  end   <- no.rain[idx + 1] - 1
  dur   <- end - start + 1
  v     <- v_all[idx]
  
  return(list("start" = start,"end" = end, "d" = dur, "v"=v))
}

aorc_rainfall_events = rain_var(data=aorc_precip_df$precip,sep_crit=6)
hist(aorc_rainfall_events$d, xlab="Event length (hrs)")
length(which(aorc_rainfall_events$d>48))/length(aorc_rainfall_events$d)
plot(aorc_rainfall_events$d,aorc_rainfall_events$v)

aorc_point_rain_sum_24hr = running_rainfall_24hr(aorc_precip_df$precip)
aorc_point_totals_24hr_window = totals_24hr_window_func(aorc_rainfall_events$start, aorc_rainfall_events$end, precip=aorc_precip_df$precip,rain_sum_24hr=aorc_point_rain_sum_24hr)
aorc_event_24hr = event_24hr_func(aorc_rainfall_events$start, aorc_rainfall_events$end, precip=aorc_precip_df$precip)


aorc_events = data.frame("start"     = aorc_rainfall_events$start, 
                         "end"       = aorc_rainfall_events$end, 
                         "dur"       = aorc_rainfall_events$d,
                         "vol"       = aorc_rainfall_events$v,
                         "sum_24"    = aorc_point_totals_24hr_window$sum_24, 
                         "sum_24_id" = aorc_point_totals_24hr_window$sum_24_id,
                         "hourly_peak"= aorc_point_totals_24hr_window$hourly_peak, 
                         "hourly_peak_id"= aorc_point_totals_24hr_window$hourly_peak_id,
                         "event_24" = aorc_event_24hr$event_24,
                         "event_24_id"= aorc_event_24hr$event_24_id)


ids <- aorc_events$event_24_id[order(aorc_events$event_24, decreasing = TRUE)[1:(5*44)]]
types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 1*24))
aorc_event_data = list(ids = ids, types = types)

ids <- aorc_events$hourly_peak_id[order(aorc_events$vol, decreasing = TRUE)[1:(5*44)]]
types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 1*24))
aorc_event_vol_data = list(ids = ids, types = types)


#Compare 24-hour totals with those from simulations
index= numeric(length(aorc_event_data$ids))
for(k in 1:length(aorc_event_data$ids)){
  index[k] = which(aorc_events$event_24_id == aorc_event_data$ids[k])
}
aorc_events$event_24[index][aorc_event_data$types=="tc"]

#No. con. on rainfall tc events expected in 44 years
length(which(df_samples$type_con_rain=="tc"))


cora.file = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")
cora_wl_df = read.csv(cora.file)[,c(1,2)]
colnames(cora_wl_df) <- c("date", "wl")

cora_wl_df$date <- as.POSIXct(round((cora_wl_df$date - 719529) * 86400 / 3600) * 3600, 
                             origin = "1970-01-01", 
                             tz = "UTC")




# ── Helpers ────────────────────────────────────────────────────────────────────
#AR2 statistic
AR2<-function(Par,Data){
  N     = length(Data)
  CDF      = sort(pgpd(q=Data,sigma=as.numeric(Par[2]),xi=as.numeric(Par[1]),u=as.numeric(Par[3])))-0.0001
  Est = N/2 + -2*sum(CDF)-sum((2-(2*(1:N)-1)/N)*log(1-CDF))
  return(Est)
}

#Generate rainfall fields
rainfall_field <- function(df_sim, rain_events, aorc_precip,
                           sim_rain_minus_ntr = NULL) {
  
  N <- nrow(df_sim)
  
  # Center a rain vector in a 145-slot window (72+offset before peak,
  # 72-offset after). Trim excess, then zero-pad to fill exactly 145.
  center_rain <- function(rain, before, after, offset = 0L) {
    left  <- 72L + offset          # target slots before peak
    right <- 72L - offset          # target slots after peak
    
    if (before > left)  rain <- rain[-(seq_len(before - left))]
    if (after  > right) rain <- rain[-seq.int(length(rain) - (after - right) + 1L,
                                              length(rain))]
    c(integer(left  - min(before, left)),
      rain,
      integer(right - min(after,  right)))
  }
  
  # ── Pre-compute all event metadata once ─────────────────────────────────────
  idx      <- match(df_sim$rain_id, rain_events$event_24_id)
  ev_start <- rain_events$start[idx]
  ev_end   <- rain_events$end[idx]
  ev_peak  <- rain_events$hourly_peak_id[idx]
  before   <- ev_peak - ev_start   # hours of rain before peak
  after    <- ev_end  - ev_peak    # hours of rain after peak
  
  is_rain_type <- df_sim$sample %in% c("r_tc", "r_non_tc")
  
  # ── Fill output matrix ───────────────────────────────────────────────────────
  sim_rain <- matrix(0L, nrow = N, ncol = 145L)
  
  for (i in seq_len(N)) {
    rain   <- df_sim$rain_factor[i] * aorc_precip[ev_start[i]:ev_end[i]]
    offset <- if (is_rain_type[i]) 0L else as.integer(sim_rain_minus_ntr[i])
    sim_rain[i, ] <- center_rain(rain, before[i], after[i], offset)
  }
  
  sim_rain
}


sim_rain = rainfall_field(df_sim=df_sim, rain_events=rain_events,aorc_precip=aorc_precip_df$precip,sim_rain_minus_ntr=sim_rain_minus_ntr)

# Weibull return-period axis for a sample of length n, scaled to 44 years
rp_axis <- function(n) (44 / n) / (1 - (1:n) / (n + 1))

# Bootstrap one event type: returns list(rp10, rp40, p_val)
bootstrap_rain <- function(obs, sim_rows, sim_rain_sums, data, perc, n_boot = 100) {
  n          <- length(obs)
  rp         <- rp_axis(n)
  rp10       <- rp40 <- p_val <- gpd_rp10 <- gpd_rp40 <- gpd_rp50 <- gpd_rp100 <- gpd_threshold <- gpd_xi <- gpd_sigma <- ar2 <-numeric(n_boot)
  sim_sorted <- vector("list", n_boot)
  p = array(0, dim=c(n_boot,n))

  for (i in seq_len(n_boot)) {
    
    samp            <- sample(length(sim_rows), n, replace = TRUE)
    s               <- sort(sim_rain_sums[sim_rows][samp])
    sim_sorted[[i]] <- s
    
    rp10[i]         <- approx(rp, s, xout = 10)$y
    rp40[i]         <- approx(rp, s, xout = 40)$y
    
    gpd = GPD_Fit(Data=s, Data_Full=data, Thres=min(s), Method="Solari")
    gpd$Rate = n/44
    
    #Extract ten year event
    gpd_rp10[i]  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10*gpd$Rate)^gpd$xi - 1)
    gpd_rp40[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((40*gpd$Rate)^gpd$xi - 1)  
    gpd_rp50[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50*gpd$Rate)^gpd$xi - 1)  
    gpd_rp100[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100*gpd$Rate)^gpd$xi - 1)  
    gpd_threshold[i]  <- gpd$Threshold 
    gpd_xi[i]  <- gpd$xi
    gpd_sigma[i] <- gpd$sigma   
    ar2[i] <- AR2(c(gpd_xi[i],gpd_sigma[i],gpd_threshold[i]),s)
    
    p_val[i]        <- ks.test(sim_rain_sums[sim_rows][samp], obs)$p.value
    
    p[i,] = quantile(sim_rain_sums[sim_rows][samp],perc)
  }
  
  p_025 = apply(p, 2, function(x) quantile(x, 0.025))
  p_05 = apply(p, 2, function(x) quantile(x, 0.05))
  p_50 = apply(p, 2, function(x) quantile(x, 0.5))
  p_95 = apply(p, 2, function(x) quantile(x, 0.95))
  p_975 = apply(p, 2, function(x) quantile(x, 0.975))
  
  list(rp10 = rp10, rp40 = rp40, 
       p_val_reject = mean(p_val < 0.05),
       gpd_rp10 = gpd_rp10, gpd_rp40 =gpd_rp40, gpd_rp50 =gpd_rp50, gpd_rp100 =gpd_rp100, 
       gpd_threshold = gpd_threshold, gpd_xi = gpd_xi, gpd_sigma =gpd_sigma, ar2 = ar2,
       p_025 = p_025, p_05 = p_05, p_50 = p_50, p_95 = p_95, p_975 = p_975,
       sim_sorted   = sim_sorted)
}


# Summarise bootstrap quantiles into a named vector
boot_quantiles <- function(x, probs = c(0.05, 0.5, 0.95),
                           suffix = "") {
  q <- quantile(x, probs)
  setNames(q, paste0("rp_", c("05","50","95"), suffix))
}

# Bootstrap Kendall correlation: returns (mean tau, rejection rate)
bootstrap_kendall <- function(sim_rows, sim_rain_sums, sim_ntr,
                              n, n_boot = 100) {
  tau  <- numeric(n_boot)
  pval <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    samp    <- sample(length(sim_rows), n, replace = TRUE)
    idx     <- sim_rows[samp]
    ct      <- cor.test(sim_rain_sums[idx], apply(sim_ntr[idx, ], 1, max),method="kendall")
    tau[i]  <- ct$estimate
    pval[i] <- ct$p.value
  }
  list(tau = mean(tau), reject_rate = mean(pval < 0.05))
}


# ── Index lookup (vectorised) ──────────────────────────────────────────────────
index <- match(aorc_event_data$ids, aorc_events$event_24_id)

# Pre-compute quantities used repeatedly
sim_rain_sums  <- rowSums(sim_rain)
is_tc          <- aorc_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")

# ── Rain marginal bootstrap ────────────────────────────────────────────────────
rain_res <- list()
for (type in c("non_tc", "tc")) {
  obs      <- aorc_events$event_24[index][if (type == "tc") is_tc else !is_tc]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- length(obs)
  rp       <- rp_axis(n)
  
  s               <- sort(obs)
  obs_rp10        <- approx(rp, s, xout = 10)$y
  obs_rp40         <- approx(rp, s, xout = 40)$y
  
  
  gpd <- GPD_Fit(Data=s, Data_Full=aorc_precip_df$precip, Thres=min(s), Method="Solari")
  gpd$Rate <- n/44
  
  #Extract ten year event
  obs_gpd_rp10  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp40  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((40*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp50  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp100 <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100*gpd$Rate)^gpd$xi - 1)  
  obs_gpd_threshold  <- gpd$Threshold 
  obs_gpd_xi  <- gpd$xi
  obs_gpd_sigma <- gpd$sigma
  obs_ar2 <- AR2(c(obs_gpd_xi,obs_gpd_sigma,obs_gpd_threshold),obs)
  
  obs_perc = ((1:n)-0.5)/n

  bt <- bootstrap_rain(obs, sim_rows, sim_rain_sums,data=aorc_precip_df$precip, perc=obs_perc)
  
  rain_res[[type]] <- list(
    obs = obs,
    obs_rp10      = obs_rp10,
    obs_rp40      = obs_rp40,
    obs_prec = obs_perc, 
    obs_gpd_rp10  = obs_gpd_rp10,
    obs_gpd_rp40  = obs_gpd_rp40,
    obs_gpd_rp50  = obs_gpd_rp50,
    obs_gpd_rp100 = obs_gpd_rp100,  

    obs_gpd_threshold  =  obs_gpd_threshold, 
    obs_gpd_xi  =  obs_gpd_xi,
    obs_gpd_sigma =  obs_gpd_sigma,   
  
    obs_ar2 =  obs_ar2,
    
    rp10_025      = quantile(bt$rp10,0.025),
    rp10_05      = quantile(bt$rp10,0.05),
    rp10_25      = quantile(bt$rp10,0.25),
    rp10_50      = quantile(bt$rp10,0.50),
    rp10_75      = quantile(bt$rp10,0.75),
    rp10_95      = quantile(bt$rp10,0.95),
    rp10_975      = quantile(bt$rp10,0.975),
    
    rp40_025      = quantile(bt$rp40,0.025),
    rp40_05      = quantile(bt$rp40,0.05),
    rp40_25      = quantile(bt$rp40,0.25),
    rp40_50      = quantile(bt$rp40,0.50),
    rp40_75      = quantile(bt$rp40,0.75),
    rp40_95      = quantile(bt$rp40,0.95),
    rp40_975      = quantile(bt$rp40,0.975),
    
    gpd_rp10_025       = quantile(bt$gpd_rp10,0.025), 
    gpd_rp10_05       = quantile(bt$gpd_rp10,0.05), 
    gpd_rp10_25       = quantile(bt$gpd_rp10,0.25), 
    gpd_rp10_50       = quantile(bt$gpd_rp10,0.50),
    gpd_rp10_75       = quantile(bt$gpd_rp10,0.75),
    gpd_rp10_95       = quantile(bt$gpd_rp10,0.95), 
    gpd_rp10_975       = quantile(bt$gpd_rp10,0.975), 
    
    gpd_rp40_025       = quantile(bt$gpd_rp40,0.025),
    gpd_rp40_05       = quantile(bt$gpd_rp40,0.05),
    gpd_rp40_25       = quantile(bt$gpd_rp40,0.25),
    gpd_rp40_50       = quantile(bt$gpd_rp40,0.50),
    gpd_rp40_75       = quantile(bt$gpd_rp40,0.75),
    gpd_rp40_95       = quantile(bt$gpd_rp40,0.95),
    gpd_rp40_975       = quantile(bt$gpd_rp40,0.975),
    
    gpd_rp50_025       = quantile(bt$gpd_rp50,0.025),
    gpd_rp50_05        = quantile(bt$gpd_rp50,0.05), 
    gpd_rp50_25       = quantile(bt$gpd_rp50,0.25), 
    gpd_rp50_50       = quantile(bt$gpd_rp50,0.50),
    gpd_rp50_75       = quantile(bt$gpd_rp50,0.75),
    gpd_rp50_95       = quantile(bt$gpd_rp50,0.95),
    gpd_rp50_975       = quantile(bt$gpd_rp50,0.975),
    
    gpd_rp100_025      = quantile(bt$gpd_rp100,0.025),
    gpd_rp100_05      = quantile(bt$gpd_rp100,0.05),
    gpd_rp100_25      = quantile(bt$gpd_rp100,0.25),
    gpd_rp100_50      = quantile(bt$gpd_rp100,0.50),
    gpd_rp100_75      = quantile(bt$gpd_rp100,0.75),
    gpd_rp100_95      = quantile(bt$gpd_rp100,0.95),
    gpd_rp100_975      = quantile(bt$gpd_rp100,0.975),
    
    gpd_threshold_025  = quantile(bt$gpd_threshold, 0.025),
    gpd_threshold_05  = quantile(bt$gpd_threshold, 0.05),
    gpd_threshold_25  = quantile(bt$gpd_threshold, 0.25),
    gpd_threshold_50  = quantile(bt$gpd_threshold, 0.50),
    gpd_threshold_75  = quantile(bt$gpd_threshold, 0.75),
    gpd_threshold_95  = quantile(bt$gpd_threshold, 0.95),
    gpd_threshold_975  = quantile(bt$gpd_threshold, 0.975),
    
    gpd_xi_025         = quantile(bt$gpd_xi,0.025),
    gpd_xi_05         = quantile(bt$gpd_xi,0.05),
    gpd_xi_25         = quantile(bt$gpd_xi,0.25),
    gpd_xi_50         = quantile(bt$gpd_xi,0.50),
    gpd_xi_75         = quantile(bt$gpd_xi,0.75),
    gpd_xi_95         = quantile(bt$gpd_xi,0.95),
    gpd_xi_975         = quantile(bt$gpd_xi,0.975),
    
    gpd_sigma_025      = quantile(bt$gpd_sigma, 0.025),
    gpd_sigma_05      = quantile(bt$gpd_sigma, 0.05),
    gpd_sigma_25      = quantile(bt$gpd_sigma, 0.25),
    gpd_sigma_50      = quantile(bt$gpd_sigma, 0.50),
    gpd_sigma_75      = quantile(bt$gpd_sigma, 0.75),
    gpd_sigma_95      = quantile(bt$gpd_sigma, 0.95),
    gpd_sigma_975      = quantile(bt$gpd_sigma, 0.975),
    
    ar2_025            = quantile(bt$ar2, 0.025),
    ar2_05            = quantile(bt$ar2, 0.05),
    ar2_25            = quantile(bt$ar2, 0.25),
    ar2_50            = quantile(bt$ar2, 0.50),
    ar2_75            = quantile(bt$ar2, 0.75),
    ar2_95            = quantile(bt$ar2, 0.95),
    ar2_975            = quantile(bt$ar2, 0.975),
    
    p_025             =  bt$p_025,
    p_05              =  bt$p_05,
    p_50              =  bt$p_50,
    p_95              =  bt$p_95,
    p_975             =  bt$p_975,
    
    p_val_reject = bt$p_val_reject
  )
}

#── Rainfall volume marginal bootstrap ────────────────────────────────────────────────────
index <- match(aorc_event_vol_data$ids, aorc_events$hourly_peak_id)

# Pre-compute quantities used repeatedly
sim_rain_sums  <- rowSums(sim_rain)
is_tc          <- aorc_event_vol_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")

rain_res_vol <- list()
for (type in c("non_tc", "tc")) {
  obs      <- aorc_events$vol[index][if (type == "tc") is_tc else !is_tc]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- length(obs)
  rp       <- rp_axis(n)
  
  s               <- sort(obs)
  obs_rp10        <- approx(rp, s, xout = 10)$y
  obs_rp40         <- approx(rp, s, xout = 40)$y
  
  
  gpd <- GPD_Fit(Data=s, Data_Full=aorc_precip_df$precip, Thres=min(s), Method="Solari")
  gpd$Rate <- n/44
  
  #Extract ten year event
  obs_gpd_rp10  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp40  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((40*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp50  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp100 <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100*gpd$Rate)^gpd$xi - 1)  
  obs_gpd_threshold  <- gpd$Threshold 
  obs_gpd_xi  <- gpd$xi
  obs_gpd_sigma <- gpd$sigma
  obs_ar2 <- AR2(c(obs_gpd_xi,obs_gpd_sigma,obs_gpd_threshold),obs)
  
  obs_perc = ((1:n)-0.5)/n
  
  bt <- bootstrap_rain(obs, sim_rows, sim_rain_sums,data=aorc_precip_df$precip, perc=obs_perc)
  
  rain_res_vol[[type]] <- list(
    obs = obs,
    obs_rp10      = obs_rp10,
    obs_rp40      = obs_rp40,
    obs_prec = obs_perc, 
    obs_gpd_rp10  = obs_gpd_rp10,
    obs_gpd_rp40  = obs_gpd_rp40,
    obs_gpd_rp50  = obs_gpd_rp50,
    obs_gpd_rp100 = obs_gpd_rp100,  
    
    obs_gpd_threshold  =  obs_gpd_threshold, 
    obs_gpd_xi  =  obs_gpd_xi,
    obs_gpd_sigma =  obs_gpd_sigma,   
    
    obs_ar2 =  obs_ar2,
    
    rp10_025      = quantile(bt$rp10,0.025),
    rp10_05      = quantile(bt$rp10,0.05),
    rp10_25      = quantile(bt$rp10,0.25),
    rp10_50      = quantile(bt$rp10,0.50),
    rp10_75      = quantile(bt$rp10,0.75),
    rp10_95      = quantile(bt$rp10,0.95),
    rp10_975      = quantile(bt$rp10,0.975),
    
    rp40_025     = quantile(bt$rp40,0.025),
    rp40_05      = quantile(bt$rp40,0.05),
    rp40_25      = quantile(bt$rp40,0.25),
    rp40_50      = quantile(bt$rp40,0.50),
    rp40_75      = quantile(bt$rp40,0.75),
    rp40_95      = quantile(bt$rp40,0.95),
    rp40_975     = quantile(bt$rp40,0.975),
    
    gpd_rp10_025      = quantile(bt$gpd_rp10,0.025), 
    gpd_rp10_05       = quantile(bt$gpd_rp10,0.05), 
    gpd_rp10_25       = quantile(bt$gpd_rp10,0.25), 
    gpd_rp10_50       = quantile(bt$gpd_rp10,0.50),
    gpd_rp10_75       = quantile(bt$gpd_rp10,0.75),
    gpd_rp10_95       = quantile(bt$gpd_rp10,0.95), 
    gpd_rp10_975      = quantile(bt$gpd_rp10,0.975), 
    
    gpd_rp40_025      = quantile(bt$gpd_rp40,0.025),
    gpd_rp40_05       = quantile(bt$gpd_rp40,0.05),
    gpd_rp40_25       = quantile(bt$gpd_rp40,0.25),
    gpd_rp40_50       = quantile(bt$gpd_rp40,0.50),
    gpd_rp40_75       = quantile(bt$gpd_rp40,0.75),
    gpd_rp40_95       = quantile(bt$gpd_rp40,0.95),
    gpd_rp40_975      = quantile(bt$gpd_rp40,0.975),
    
    gpd_rp50_025      = quantile(bt$gpd_rp50,0.025),
    gpd_rp50_05       = quantile(bt$gpd_rp50,0.05), 
    gpd_rp50_25       = quantile(bt$gpd_rp50,0.25), 
    gpd_rp50_50       = quantile(bt$gpd_rp50,0.50),
    gpd_rp50_75       = quantile(bt$gpd_rp50,0.75),
    gpd_rp50_95       = quantile(bt$gpd_rp50,0.95),
    gpd_rp50_975      = quantile(bt$gpd_rp50,0.975),
    
    gpd_rp100_025     = quantile(bt$gpd_rp100,0.025),
    gpd_rp100_05      = quantile(bt$gpd_rp100,0.05),
    gpd_rp100_25      = quantile(bt$gpd_rp100,0.25),
    gpd_rp100_50      = quantile(bt$gpd_rp100,0.50),
    gpd_rp100_75      = quantile(bt$gpd_rp100,0.75),
    gpd_rp100_95      = quantile(bt$gpd_rp100,0.95),
    gpd_rp100_975     = quantile(bt$gpd_rp100,0.975),
    
    gpd_threshold_025  = quantile(bt$gpd_threshold, 0.025),
    gpd_threshold_05  = quantile(bt$gpd_threshold, 0.05),
    gpd_threshold_25  = quantile(bt$gpd_threshold, 0.25),
    gpd_threshold_50  = quantile(bt$gpd_threshold, 0.50),
    gpd_threshold_75  = quantile(bt$gpd_threshold, 0.75),
    gpd_threshold_95  = quantile(bt$gpd_threshold, 0.95),
    gpd_threshold_975 = quantile(bt$gpd_threshold, 0.975),
    
    gpd_xi_025        = quantile(bt$gpd_xi,0.025),
    gpd_xi_05         = quantile(bt$gpd_xi,0.05),
    gpd_xi_25         = quantile(bt$gpd_xi,0.25),
    gpd_xi_50         = quantile(bt$gpd_xi,0.50),
    gpd_xi_75         = quantile(bt$gpd_xi,0.75),
    gpd_xi_95         = quantile(bt$gpd_xi,0.95),
    gpd_xi_975        = quantile(bt$gpd_xi,0.975),
    
    gpd_sigma_025     = quantile(bt$gpd_sigma, 0.025),
    gpd_sigma_05      = quantile(bt$gpd_sigma, 0.05),
    gpd_sigma_25      = quantile(bt$gpd_sigma, 0.25),
    gpd_sigma_50      = quantile(bt$gpd_sigma, 0.50),
    gpd_sigma_75      = quantile(bt$gpd_sigma, 0.75),
    gpd_sigma_95      = quantile(bt$gpd_sigma, 0.95),
    gpd_sigma_975     = quantile(bt$gpd_sigma, 0.975),
    
    ar2_025           = quantile(bt$ar2, 0.025),
    ar2_05            = quantile(bt$ar2, 0.05),
    ar2_25            = quantile(bt$ar2, 0.25),
    ar2_50            = quantile(bt$ar2, 0.50),
    ar2_75            = quantile(bt$ar2, 0.75),
    ar2_95            = quantile(bt$ar2, 0.95),
    ar2_975           = quantile(bt$ar2, 0.975),
    
    p_025            =  bt$p_025,
    p_05             =  bt$p_05,
    p_50             =  bt$p_50,
    p_95             =  bt$p_95,
    p_975            =  bt$p_975,
    
    p_val_reject = bt$p_val_reject
  )
}


# ── Kendall correlation bootstrap con rainfall ─────────────────────────────────────────────
library(energy)

sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
is_tc          <- aorc_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")


cor_con_rain_res = list()
for (type in c("tc", "non_tc")) {

  obs_df      <- df_samples[df_samples$type_con_rain==type,c('rain_con_rain','ntr_con_rain')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_sums, wl$sim_ntr_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_sums, sim_ntr_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_con_rain_res[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
}



# ── Kendall correlation bootstrap  con ntr─────────────────────────────────────────────
library(energy)

sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")


cor_con_ntr_res = list()
for (type in c("tc", "non_tc")) {
  obs_df      <- df_samples[df_samples$type_con_ntr==type,c('rain_con_ntr','ntr_con_ntr')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc

  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_sums, wl$sim_ntr_trim,
                                   n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_sums, sim_ntr_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_con_ntr_res[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
  
}
res_rain_list[[i]] <- data.frame(
  obs_rp10_tc           = rain_res$tc$obs_rp10,
  obs_rp40_tc           = rain_res$tc$obs_rp40,
  obs_gpd_rp10_tc       = rain_res$tc$obs_gpd_rp10,
  obs_gpd_rp40_tc       = rain_res$tc$obs_gpd_rp40,
  obs_gpd_rp50_tc       = rain_res$tc$obs_gpd_rp50,
  obs_gpd_rp100_tc      = rain_res$tc$obs_gpd_rp100,
  
  obs_gpd_threshold_tc     = rain_res$tc$obs_gpd_threshold, 
  obs_gpd_xi_tc            = rain_res$tc$obs_gpd_xi,
  obs_gpd_sigma_tc         = rain_res$tc$obs_gpd_sigma,   
  obs_ar2_tc               = rain_res$tc$obs_ar2,

  rp10_025_tc       = as.numeric(rain_res$tc$rp10_025),
  rp10_05_tc       = as.numeric(rain_res$tc$rp10_05),
  rp10_25_tc       = as.numeric(rain_res$tc$rp10_25),
  rp10_50_tc       = as.numeric(rain_res$tc$rp10_50),
  rp10_75_tc       = as.numeric(rain_res$tc$rp10_75),
  rp10_95_tc       = as.numeric(rain_res$tc$rp10_95),
  rp10_975_tc       = as.numeric(rain_res$tc$rp10_975),
  
  rp40_025_tc       = rain_res$tc$rp40_025,
  rp40_05_tc       = rain_res$tc$rp40_05,
  rp40_25_tc       = rain_res$tc$rp40_25,
  rp40_50_tc       = rain_res$tc$rp40_50,
  rp40_75_tc       = rain_res$tc$rp40_75,
  rp40_95_tc       = rain_res$tc$rp40_95,
  rp40_975_tc       = rain_res$tc$rp40_975,
  
  gpd_rp10_025_tc       = rain_res$tc$gpd_rp10_025,
  gpd_rp10_05_tc       = rain_res$tc$gpd_rp10_05,
  gpd_rp10_25_tc       = rain_res$tc$gpd_rp10_25,
  gpd_rp10_50_tc       = rain_res$tc$gpd_rp10_50,
  gpd_rp10_75_tc       = rain_res$tc$gpd_rp10_75,
  gpd_rp10_95_tc       = rain_res$tc$gpd_rp10_95,
  gpd_rp10_975_tc       = rain_res$tc$gpd_rp10_975,
  
  gpd_rp40_025_tc       = rain_res$tc$gpd_rp40_025,
  gpd_rp40_05_tc       = rain_res$tc$gpd_rp40_05,
  gpd_rp40_25_tc       = rain_res$tc$gpd_rp40_25,
  gpd_rp40_50_tc       = rain_res$tc$gpd_rp40_50,
  gpd_rp40_75_tc       = rain_res$tc$gpd_rp40_75,
  gpd_rp40_95_tc       = rain_res$tc$gpd_rp40_95,
  gpd_rp40_975_tc       = rain_res$tc$gpd_rp40_975,
  
  gpd_rp50_025_tc       = rain_res$tc$gpd_rp50_025,
  gpd_rp50_05_tc       = rain_res$tc$gpd_rp50_05,
  gpd_rp50_25_tc       = rain_res$tc$gpd_rp50_25,
  gpd_rp50_50_tc       = rain_res$tc$gpd_rp50_50,
  gpd_rp50_75_tc       = rain_res$tc$gpd_rp50_75,
  gpd_rp50_95_tc       = rain_res$tc$gpd_rp50_95,
  gpd_rp50_975_tc       = rain_res$tc$gpd_rp50_975,
  
  gpd_rp100_025_tc      = rain_res$tc$gpd_rp100_025,
  gpd_rp100_05_tc       = rain_res$tc$gpd_rp100_05,
  gpd_rp100_25_tc       = rain_res$tc$gpd_rp100_25,
  gpd_rp100_50_tc       = rain_res$tc$gpd_rp100_50,
  gpd_rp100_75_tc       = rain_res$tc$gpd_rp100_75,
  gpd_rp100_95_tc       = rain_res$tc$gpd_rp100_95,
  gpd_rp100_975_tc      = rain_res$tc$gpd_rp100_975,
  
  gpd_threshold_025_tc           = rain_res$tc$gpd_threshold_025,
  gpd_threshold_05_tc            = rain_res$tc$gpd_threshold_05,
  gpd_threshold_25_tc            = rain_res$tc$gpd_threshold_25,
  gpd_threshold_50_tc            = rain_res$tc$gpd_threshold_50,
  gpd_threshold_75_tc            = rain_res$tc$gpd_threshold_75,
  gpd_threshold_95_tc            = rain_res$tc$gpd_threshold_95,
  gpd_threshold_975_tc           = rain_res$tc$gpd_threshold_975,
  
  gpd_xi_025_tc                  = rain_res$tc$gpd_xi_025,
  gpd_xi_05_tc                   = rain_res$tc$gpd_xi_05,
  gpd_xi_25_tc                   = rain_res$tc$gpd_xi_25,
  gpd_xi_50_tc                   = rain_res$tc$gpd_xi_50,
  gpd_xi_75_tc                   = rain_res$tc$gpd_xi_75,
  gpd_xi_95_tc                   = rain_res$tc$gpd_xi_95,
  gpd_xi_975_tc                  = rain_res$tc$gpd_xi_975,
  
  gpd_sigma_025_tc               = rain_res$tc$gpd_sigma_025,
  gpd_sigma_05_tc                = rain_res$tc$gpd_sigma_05,
  gpd_sigma_25_tc                = rain_res$tc$gpd_sigma_25,
  gpd_sigma_50_tc                = rain_res$tc$gpd_sigma_50,
  gpd_sigma_75_tc                = rain_res$tc$gpd_sigma_75,
  gpd_sigma_95_tc                = rain_res$tc$gpd_sigma_95,
  gpd_sigma_975_tc               = rain_res$tc$gpd_sigma_975,
  
  ar2_025_tc                     = rain_res$tc$ar2_025,
  ar2_05_tc                      = rain_res$tc$ar2_05,
  ar2_25_tc                      = rain_res$tc$ar2_25,
  ar2_50_tc                      = rain_res$tc$ar2_50,
  ar2_75_tc                      = rain_res$tc$ar2_75,
  ar2_95_tc                      = rain_res$tc$ar2_95,
  ar2_975_tc                     = rain_res$tc$ar2_975,
  
  p_val_reject_tc             = rain_res$tc$p_val_reject,
  
  obs_rp10_non_tc           = rain_res$non_tc$obs_rp10,
  obs_rp40_non_tc           = rain_res$non_tc$obs_rp40,
  obs_gpd_rp10_non_tc       = rain_res$non_tc$obs_gpd_rp10,
  obs_gpd_rp40_non_tc       = rain_res$non_tc$obs_gpd_rp40,
  obs_gpd_rp50_non_tc       = rain_res$non_tc$obs_gpd_rp50,
  obs_gpd_rp100_non_tc      = rain_res$non_tc$obs_gpd_rp100,
  
  obs_gpd_threshold_non_tc     = rain_res$non_tc$obs_gpd_threshold, 
  obs_gpd_xi_non_tc            = rain_res$non_tc$obs_gpd_xi,
  obs_gpd_sigma_non_tc         = rain_res$non_tc$obs_gpd_sigma,   
  obs_ar2_non_tc               = rain_res$non_tc$obs_ar2,
  
  rp10_025_non_tc      = rain_res$non_tc$rp10_025,
  rp10_05_non_tc       = rain_res$non_tc$rp10_05,
  rp10_25_non_tc       = rain_res$non_tc$rp10_25,
  rp10_50_non_tc       = rain_res$non_tc$rp10_50,
  rp10_75_non_tc       = rain_res$non_tc$rp10_75,
  rp10_95_non_tc       = rain_res$non_tc$rp10_95,
  rp10_975_non_tc      = rain_res$non_tc$rp10_975,
  
  rp40_025_non_tc      = rain_res$non_tc$rp40_025,
  rp40_05_non_tc       = rain_res$non_tc$rp40_05,
  rp40_25_non_tc       = rain_res$non_tc$rp40_25,
  rp40_50_non_tc       = rain_res$non_tc$rp40_50,
  rp40_75_non_tc       = rain_res$non_tc$rp40_75,
  rp40_95_non_tc       = rain_res$non_tc$rp40_95,
  rp40_975_non_tc      = rain_res$non_tc$rp40_975,
  
  gpd_rp10_025_non_tc       = rain_res$non_tc$gpd_rp10_025,
  gpd_rp10_05_non_tc       = rain_res$non_tc$gpd_rp10_05,
  gpd_rp10_25_non_tc       = rain_res$non_tc$gpd_rp10_25,
  gpd_rp10_50_non_tc       = rain_res$non_tc$gpd_rp10_50,
  gpd_rp10_75_non_tc       = rain_res$non_tc$gpd_rp10_75,
  gpd_rp10_95_non_tc       = rain_res$non_tc$gpd_rp10_95,
  gpd_rp10_975_non_tc       = rain_res$non_tc$gpd_rp10_975,
  
  gpd_rp40_025_non_tc       = rain_res$non_tc$gpd_rp40_025,
  gpd_rp40_05_non_tc       = rain_res$non_tc$gpd_rp40_05,
  gpd_rp40_25_non_tc       = rain_res$non_tc$gpd_rp40_25,
  gpd_rp40_50_non_tc       = rain_res$non_tc$gpd_rp40_50,
  gpd_rp40_75_non_tc       = rain_res$non_tc$gpd_rp40_75,
  gpd_rp40_95_non_tc       = rain_res$non_tc$gpd_rp40_95,
  gpd_rp40_975_non_tc       = rain_res$non_tc$gpd_rp40_975,
  
  gpd_rp50_025_non_tc       = rain_res$non_tc$gpd_rp50_025,
  gpd_rp50_05_non_tc       = rain_res$non_tc$gpd_rp50_05,
  gpd_rp50_25_non_tc       = rain_res$non_tc$gpd_rp50_25,
  gpd_rp50_50_non_tc       = rain_res$non_tc$gpd_rp50_50,
  gpd_rp50_75_non_tc       = rain_res$non_tc$gpd_rp50_75,
  gpd_rp50_95_non_tc       = rain_res$non_tc$gpd_rp50_95,
  gpd_rp50_975_non_tc       = rain_res$non_tc$gpd_rp50_975,

  gpd_rp100_025_non_tc      = rain_res$non_tc$gpd_rp100_025,
  gpd_rp100_05_non_tc      = rain_res$non_tc$gpd_rp100_05,
  gpd_rp100_25_non_tc      = rain_res$non_tc$gpd_rp100_25,
  gpd_rp100_50_non_tc      = rain_res$non_tc$gpd_rp100_50,
  gpd_rp100_75_non_tc      = rain_res$non_tc$gpd_rp100_75,
  gpd_rp100_95_non_tc      = rain_res$non_tc$gpd_rp100_95,
  gpd_rp100_975_non_tc      = rain_res$non_tc$gpd_rp100_975,
  
  gpd_threshold_025_non_tc            = rain_res$non_tc$gpd_threshold_025,
  gpd_threshold_05_non_tc            = rain_res$non_tc$gpd_threshold_05,
  gpd_threshold_25_non_tc            = rain_res$non_tc$gpd_threshold_25,
  gpd_threshold_50_non_tc            = rain_res$non_tc$gpd_threshold_50,
  gpd_threshold_75_non_tc            = rain_res$non_tc$gpd_threshold_75,
  gpd_threshold_95_non_tc            = rain_res$non_tc$gpd_threshold_95,
  gpd_threshold_975_non_tc            = rain_res$non_tc$gpd_threshold_975,
  
  gpd_xi_025_non_tc                   = rain_res$non_tc$gpd_xi_025,
  gpd_xi_05_non_tc                   = rain_res$non_tc$gpd_xi_05,
  gpd_xi_25_non_tc                   = rain_res$non_tc$gpd_xi_25,
  gpd_xi_50_non_tc                   = rain_res$non_tc$gpd_xi_50,
  gpd_xi_75_non_tc                   = rain_res$non_tc$gpd_xi_75,
  gpd_xi_95_non_tc                   = rain_res$non_tc$gpd_xi_95,
  gpd_xi_975_non_tc                   = rain_res$non_tc$gpd_xi_975,
  
  gpd_sigma_025_non_tc                = rain_res$non_tc$gpd_sigma_025,
  gpd_sigma_05_non_tc                = rain_res$non_tc$gpd_sigma_05,
  gpd_sigma_25_non_tc                = rain_res$non_tc$gpd_sigma_25,
  gpd_sigma_50_non_tc                = rain_res$non_tc$gpd_sigma_50,
  gpd_sigma_75_non_tc                = rain_res$non_tc$gpd_sigma_75,
  gpd_sigma_95_non_tc                = rain_res$non_tc$gpd_sigma_95,
  gpd_sigma_975_non_tc                = rain_res$non_tc$gpd_sigma_975,
  
  ar2_025_non_tc                      = rain_res$non_tc$ar2_025,
  ar2_05_non_tc                      = rain_res$non_tc$ar2_05,
  ar2_25_non_tc                      = rain_res$non_tc$ar2_25,
  ar2_50_non_tc                      = rain_res$non_tc$ar2_50,
  ar2_75_non_tc                      = rain_res$non_tc$ar2_75,
  ar2_95_non_tc                      = rain_res$non_tc$ar2_95,
  ar2_975_non_tc                      = rain_res$non_tc$ar2_975,
  
  p_val_reject_non_tc             = rain_res$non_tc$p_val_reject
)


res_rain_df_tc_list[[i]] = data.frame(aorc = rep(i,length(rain_res$tc$p_05)),
                                      obs_tc            =  rain_res$tc$obs,
                                      p025_tc            =  rain_res$tc$p_025,
                                      p05_tc            =  rain_res$tc$p_05,
                                      p50_tc            =  rain_res$tc$p_50,
                                      p95_tc            =  rain_res$tc$p_95,
                                      p975_tc            =  rain_res$tc$p_975)

res_rain_df_non_tc_list[[i]] = data.frame(aorc = rep(i,length(rain_res$non_tc$p_05)),
                                          obs_non_tc            =  rain_res$non_tc$obs,
                                          p025_non_tc            =  rain_res$non_tc$p_025,
                                          p05_non_tc            =  rain_res$non_tc$p_05,
                                          p50_non_tc            =  rain_res$non_tc$p_50,
                                          p95_non_tc            =  rain_res$non_tc$p_95,
                                          p975_non_tc            =  rain_res$non_tc$p_975)



res_vol_list[[i]] <- data.frame(
  obs_rp10_tc           = rain_res_vol$tc$obs_rp10,
  obs_rp40_tc           = rain_res_vol$tc$obs_rp40,
  obs_gpd_rp10_tc       = rain_res_vol$tc$obs_gpd_rp10,
  obs_gpd_rp40_tc       = rain_res_vol$tc$obs_gpd_rp40,
  obs_gpd_rp50_tc       = rain_res_vol$tc$obs_gpd_rp50,
  obs_gpd_rp100_tc      = rain_res_vol$tc$obs_gpd_rp100,
  
  obs_gpd_threshold_tc     = rain_res_vol$tc$obs_gpd_threshold, 
  obs_gpd_xi_tc            = rain_res_vol$tc$obs_gpd_xi,
  obs_gpd_sigma_tc         = rain_res_vol$tc$obs_gpd_sigma,   
  obs_ar2_tc               = rain_res_vol$tc$obs_ar2,
  
  rp10_025_tc       = rain_res_vol$tc$rp10_025,
  rp10_05_tc       = rain_res_vol$tc$rp10_05,
  rp10_25_tc       = rain_res_vol$tc$rp10_25,
  rp10_50_tc       = rain_res_vol$tc$rp10_50,
  rp10_75_tc       = rain_res_vol$tc$rp10_75,
  rp10_95_tc       = rain_res_vol$tc$rp10_95,
  rp10_975_tc       = rain_res_vol$tc$rp10_975,
  
  rp40_025_tc       = rain_res_vol$tc$rp40_025,
  rp40_05_tc       = rain_res_vol$tc$rp40_05,
  rp40_25_tc       = rain_res_vol$tc$rp40_25,
  rp40_50_tc       = rain_res_vol$tc$rp40_50,
  rp40_75_tc       = rain_res_vol$tc$rp40_75,
  rp40_95_tc       = rain_res_vol$tc$rp40_95,
  rp40_975_tc       = rain_res_vol$tc$rp40_975,

  gpd_rp10_025_tc       = rain_res_vol$tc$gpd_rp10_025,  
  gpd_rp10_05_tc       = rain_res_vol$tc$gpd_rp10_05,
  gpd_rp10_25_tc       = rain_res_vol$tc$gpd_rp10_25,
  gpd_rp10_50_tc       = rain_res_vol$tc$gpd_rp10_50,
  gpd_rp10_75_tc       = rain_res_vol$tc$gpd_rp10_75,
  gpd_rp10_95_tc       = rain_res_vol$tc$gpd_rp10_95,
  gpd_rp10_975_tc       = rain_res_vol$tc$gpd_rp10_975,
  
  gpd_rp40_025_tc       = rain_res_vol$tc$gpd_rp40_025,
  gpd_rp40_05_tc       = rain_res_vol$tc$gpd_rp40_05,
  gpd_rp40_25_tc       = rain_res_vol$tc$gpd_rp40_25,
  gpd_rp40_50_tc       = rain_res_vol$tc$gpd_rp40_50,
  gpd_rp40_75_tc       = rain_res_vol$tc$gpd_rp40_75,
  gpd_rp40_95_tc       = rain_res_vol$tc$gpd_rp40_95,
  gpd_rp40_975_tc       = rain_res_vol$tc$gpd_rp40_975,
  
  gpd_rp50_025_tc       = rain_res_vol$tc$gpd_rp50_025,
  gpd_rp50_05_tc       = rain_res_vol$tc$gpd_rp50_05,
  gpd_rp50_25_tc       = rain_res_vol$tc$gpd_rp50_25,
  gpd_rp50_50_tc       = rain_res_vol$tc$gpd_rp50_50,
  gpd_rp50_75_tc       = rain_res_vol$tc$gpd_rp50_75,
  gpd_rp50_95_tc       = rain_res_vol$tc$gpd_rp50_95,
  gpd_rp50_975_tc       = rain_res_vol$tc$gpd_rp50_975,
  
  gpd_rp100_025_tc      = rain_res_vol$tc$gpd_rp100_025,
  gpd_rp100_05_tc      = rain_res_vol$tc$gpd_rp100_05,
  gpd_rp100_25_tc      = rain_res_vol$tc$gpd_rp100_25,
  gpd_rp100_50_tc      = rain_res_vol$tc$gpd_rp100_50,
  gpd_rp100_75_tc      = rain_res_vol$tc$gpd_rp100_75,
  gpd_rp100_95_tc      = rain_res_vol$tc$gpd_rp100_95,
  gpd_rp100_975_tc      = rain_res_vol$tc$gpd_rp100_975,
  
  gpd_threshold_025_tc            = rain_res_vol$tc$gpd_threshold_025,
  gpd_threshold_05_tc            = rain_res_vol$tc$gpd_threshold_05,
  gpd_threshold_25_tc            = rain_res_vol$tc$gpd_threshold_25,
  gpd_threshold_50_tc            = rain_res_vol$tc$gpd_threshold_50,
  gpd_threshold_75_tc            = rain_res_vol$tc$gpd_threshold_75,
  gpd_threshold_95_tc            = rain_res_vol$tc$gpd_threshold_95,
  gpd_threshold_975_tc            = rain_res_vol$tc$gpd_threshold_975,
  
  gpd_xi_025_tc                   = rain_res_vol$tc$gpd_xi_025,
  gpd_xi_05_tc                   = rain_res_vol$tc$gpd_xi_05,
  gpd_xi_25_tc                   = rain_res_vol$tc$gpd_xi_25,
  gpd_xi_50_tc                   = rain_res_vol$tc$gpd_xi_50,
  gpd_xi_75_tc                   = rain_res_vol$tc$gpd_xi_75,
  gpd_xi_95_tc                   = rain_res_vol$tc$gpd_xi_95,
  gpd_xi_975_tc                   = rain_res_vol$tc$gpd_xi_975,
  
  gpd_sigma_025_tc                = rain_res_vol$tc$gpd_sigma_025,
  gpd_sigma_05_tc                = rain_res_vol$tc$gpd_sigma_05,
  gpd_sigma_25_tc                = rain_res_vol$tc$gpd_sigma_25,
  gpd_sigma_50_tc                = rain_res_vol$tc$gpd_sigma_50,
  gpd_sigma_75_tc                = rain_res_vol$tc$gpd_sigma_75,
  gpd_sigma_95_tc                = rain_res_vol$tc$gpd_sigma_95,
  gpd_sigma_975_tc                = rain_res_vol$tc$gpd_sigma_975,
  
  ar2_025_tc                      = rain_res_vol$tc$ar2_025,
  ar2_05_tc                      = rain_res_vol$tc$ar2_05,
  ar2_25_tc                      = rain_res_vol$tc$ar2_25,
  ar2_50_tc                      = rain_res_vol$tc$ar2_50,
  ar2_75_tc                      = rain_res_vol$tc$ar2_75,
  ar2_95_tc                      = rain_res_vol$tc$ar2_95,
  ar2_975_tc                      = rain_res_vol$tc$ar2_975,

  p_val_reject_tc             = rain_res_vol$tc$p_val_reject,
  
  obs_rp10_non_tc           = rain_res_vol$non_tc$obs_rp10,
  obs_rp40_non_tc           = rain_res_vol$non_tc$obs_rp40,
  obs_gpd_rp10_non_tc       = rain_res_vol$non_tc$obs_gpd_rp10,
  obs_gpd_rp40_non_tc       = rain_res_vol$non_tc$obs_gpd_rp40,
  obs_gpd_rp50_non_tc       = rain_res_vol$non_tc$obs_gpd_rp50,
  obs_gpd_rp100_non_tc      = rain_res_vol$non_tc$obs_gpd_rp100,
  
  obs_gpd_threshold_non_tc     = rain_res_vol$non_tc$obs_gpd_threshold, 
  obs_gpd_xi_non_tc            = rain_res_vol$non_tc$obs_gpd_xi,
  obs_gpd_sigma_non_tc         = rain_res_vol$non_tc$obs_gpd_sigma,   
  obs_ar2_non_tc               = rain_res_vol$non_tc$obs_ar2,
  
  
  rp10_025_non_tc       = rain_res_vol$non_tc$rp10_025,
  rp10_05_non_tc       = rain_res_vol$non_tc$rp10_05,
  rp10_25_non_tc       = rain_res_vol$non_tc$rp10_25,
  rp10_50_non_tc       = rain_res_vol$non_tc$rp10_50,
  rp10_75_non_tc       = rain_res_vol$non_tc$rp10_75,
  rp10_95_non_tc       = rain_res_vol$non_tc$rp10_95,
  rp10_975_non_tc       = rain_res_vol$non_tc$rp10_975,
  
  rp40_025_non_tc       = rain_res_vol$non_tc$rp40_025,
  rp40_05_non_tc       = rain_res_vol$non_tc$rp40_05,
  rp40_25_non_tc       = rain_res_vol$non_tc$rp40_25,
  rp40_50_non_tc       = rain_res_vol$non_tc$rp40_50,
  rp40_75_non_tc       = rain_res_vol$non_tc$rp40_75,
  rp40_95_non_tc       = rain_res_vol$non_tc$rp40_95,
  rp40_975_non_tc       = rain_res_vol$non_tc$rp40_975,
  
  gpd_rp10_025_non_tc       = rain_res_vol$non_tc$gpd_rp10_025,
  gpd_rp10_05_non_tc       = rain_res_vol$non_tc$gpd_rp10_05,
  gpd_rp10_25_non_tc       = rain_res_vol$non_tc$gpd_rp10_25,
  gpd_rp10_50_non_tc       = rain_res_vol$non_tc$gpd_rp10_50,
  gpd_rp10_75_non_tc       = rain_res_vol$non_tc$gpd_rp10_75,
  gpd_rp10_95_non_tc       = rain_res_vol$non_tc$gpd_rp10_95,
  gpd_rp10_975_non_tc       = rain_res_vol$non_tc$gpd_rp10_975,
  
  gpd_rp40_025_non_tc       = rain_res_vol$non_tc$gpd_rp40_025,
  gpd_rp40_05_non_tc       = rain_res_vol$non_tc$gpd_rp40_05,
  gpd_rp40_25_non_tc       = rain_res_vol$non_tc$gpd_rp40_25,
  gpd_rp40_50_non_tc       = rain_res_vol$non_tc$gpd_rp40_50,
  gpd_rp40_75_non_tc       = rain_res_vol$non_tc$gpd_rp40_75,
  gpd_rp40_95_non_tc       = rain_res_vol$non_tc$gpd_rp40_95,
  gpd_rp40_975_non_tc       = rain_res_vol$non_tc$gpd_rp40_975,
  
  gpd_rp50_025_non_tc       = rain_res_vol$non_tc$gpd_rp50_025,
  gpd_rp50_05_non_tc       = rain_res_vol$non_tc$gpd_rp50_05,
  gpd_rp50_25_non_tc       = rain_res_vol$non_tc$gpd_rp50_25,
  gpd_rp50_50_non_tc       = rain_res_vol$non_tc$gpd_rp50_50,
  gpd_rp50_75_non_tc       = rain_res_vol$non_tc$gpd_rp50_75,
  gpd_rp50_95_non_tc       = rain_res_vol$non_tc$gpd_rp50_95,
  gpd_rp50_975_non_tc       = rain_res_vol$non_tc$gpd_rp50_975,
  
  gpd_rp100_025_non_tc      = rain_res_vol$non_tc$gpd_rp100_025,
  gpd_rp100_05_non_tc      = rain_res_vol$non_tc$gpd_rp100_05,
  gpd_rp100_25_non_tc      = rain_res_vol$non_tc$gpd_rp100_25,
  gpd_rp100_50_non_tc      = rain_res_vol$non_tc$gpd_rp100_50,
  gpd_rp100_75_non_tc      = rain_res_vol$non_tc$gpd_rp100_75,
  gpd_rp100_95_non_tc      = rain_res_vol$non_tc$gpd_rp100_95,
  gpd_rp100_975_non_tc      = rain_res_vol$non_tc$gpd_rp100_975,
  
  gpd_threshold_025_non_tc            = rain_res_vol$non_tc$gpd_threshold_025,
  gpd_threshold_05_non_tc            = rain_res_vol$non_tc$gpd_threshold_05,
  gpd_threshold_25_non_tc            = rain_res_vol$non_tc$gpd_threshold_25,
  gpd_threshold_50_non_tc            = rain_res_vol$non_tc$gpd_threshold_50,
  gpd_threshold_75_non_tc            = rain_res_vol$non_tc$gpd_threshold_75,
  gpd_threshold_95_non_tc            = rain_res_vol$non_tc$gpd_threshold_95,
  gpd_threshold_975_non_tc            = rain_res_vol$non_tc$gpd_threshold_975,
  
  gpd_xi_025_non_tc                   = rain_res_vol$non_tc$gpd_xi_025,
  gpd_xi_05_non_tc                   = rain_res_vol$non_tc$gpd_xi_05,
  gpd_xi_25_non_tc                   = rain_res_vol$non_tc$gpd_xi_25,
  gpd_xi_50_non_tc                   = rain_res_vol$non_tc$gpd_xi_50,
  gpd_xi_75_non_tc                   = rain_res_vol$non_tc$gpd_xi_75,
  gpd_xi_95_non_tc                   = rain_res_vol$non_tc$gpd_xi_95,
  gpd_xi_975_non_tc                   = rain_res_vol$non_tc$gpd_xi_975,
  
  gpd_sigma_025_non_tc                = rain_res_vol$non_tc$gpd_sigma_025,
  gpd_sigma_05_non_tc                = rain_res_vol$non_tc$gpd_sigma_05,
  gpd_sigma_25_non_tc                = rain_res_vol$non_tc$gpd_sigma_25,
  gpd_sigma_50_non_tc                = rain_res_vol$non_tc$gpd_sigma_50,
  gpd_sigma_75_non_tc                = rain_res_vol$non_tc$gpd_sigma_75,
  gpd_sigma_95_non_tc                = rain_res_vol$non_tc$gpd_sigma_95,
  gpd_sigma_975_non_tc                = rain_res_vol$non_tc$gpd_sigma_975,
  
  ar2_025_non_tc                      = rain_res_vol$non_tc$ar2_025,
  ar2_05_non_tc                      = rain_res_vol$non_tc$ar2_05,
  ar2_25_non_tc                      = rain_res_vol$non_tc$ar2_25,
  ar2_50_non_tc                      = rain_res_vol$non_tc$ar2_50,
  ar2_75_non_tc                      = rain_res_vol$non_tc$ar2_75,
  ar2_95_non_tc                      = rain_res_vol$non_tc$ar2_95,
  ar2_975_non_tc                      = rain_res_vol$non_tc$ar2_975,

  p_val_reject_non_tc             = rain_res_vol$non_tc$p_val_reject
)


res_vol_df_tc_list[[i]] = data.frame(aorc              =  rep(i,length(rain_res_vol$tc$p_05)),
                                     obs_tc            =  rain_res_vol$tc$obs,
                                     p025_tc            =  rain_res_vol$tc$p_025,
                                     p05_tc            =  rain_res_vol$tc$p_05,
                                     p50_tc            =  rain_res_vol$tc$p_50,
                                     p95_tc            =  rain_res_vol$tc$p_95,
                                     p975_tc            =  rain_res_vol$tc$p_975)

res_vol_df_non_tc_list[[i]] = data.frame(aorc                  =  rep(i,length(rain_res_vol$non_tc$p_05)),
                                         obs_non_tc            =  rain_res_vol$non_tc$obs,
                                         p025_non_tc            =  rain_res_vol$non_tc$p_025,
                                         p05_non_tc            =  rain_res_vol$non_tc$p_05,
                                         p50_non_tc            =  rain_res_vol$non_tc$p_50,
                                         p95_non_tc            =  rain_res_vol$non_tc$p_95,
                                         p975_non_tc            =  rain_res_vol$non_tc$p_975)



res_cor_con_rain_list[[i]] <- data.frame(
  cor_r_tc = cor_con_rain_res$tc$tau,
  cor_r_tc_reject_rate = cor_con_rain_res$tc$reject_rate,
  cor_r_non_tc = cor_con_rain_res$non_tc$tau,
  cor_r_non_tc_reject_rate = cor_con_rain_res$non_tc$reject_rate)

#cor_con_rain_tc_edist_test_val = cor_con_rain_res$tc$edist_test_val,
#cor_con_rain_tc_edist_test_p_val = cor_con_rain_res$tc$edist_test_p_val,
#cor_con_rain_non_tc_edist_test_val = cor_con_rain_res$non_tc$edist_test_val,
#cor_con_rain_non_tc_edist_test_p_val = cor_con_rain_res$non_tc$edist_test_p_val

res_cor_con_ntr_list[[i]] <- data.frame(
  cor_ntr_tc = cor_con_ntr_res$tc$tau,
  cor_ntr_tc_reject_rate = cor_con_ntr_res$tc$reject_rate,
  cor_ntr_non_tc = cor_con_ntr_res$non_tc$tau,
  cor_ntr_non_tc_reject_rate = cor_con_ntr_res$non_tc$reject_rate)



# ── Kendall correlation bootstrap con rainfall ─────────────────────────────────────────────
library(energy)

sim_rain_sums  <- rowSums(sim_rain)
sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")


cor_rain_wl_con_rain_res = list()
for (type in c("tc", "non_tc")) {
  
  obs_df      <- df_samples[df_samples$type_con_rain==type,c('rain_con_rain','ntr_con_rain')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_sums, wl$sim_wl_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_sums, sim_wl_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_rain_wl_con_rain_res[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
}

res_cor_rain_wl_con_rain_list[[i]] <- data.frame(
  cor_r_tc = cor_rain_wl_con_rain_res$tc$tau,
  cor_r_tc_reject_rate = cor_rain_wl_con_rain_res$tc$reject_rate,
  cor_r_non_tc = cor_rain_wl_con_rain_res$non_tc$tau,
  cor_r_non_tc_reject_rate = cor_rain_wl_con_rain_res$non_tc$reject_rate)


sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")


cor_rain_wl_con_wl_res = list()
for (type in c("tc", "non_tc")) {
  obs_df      <- df_samples[df_samples$type_con_ntr==type,c('rain_con_ntr','ntr_con_ntr')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_sums, wl$sim_wl_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_sums, sim_wl_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_rain_wl_con_wl_res[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
}

res_cor_rain_wl_con_wl_list[[i]] <- data.frame(
  cor_r_tc = cor_rain_wl_con_wl_res$tc$tau,
  cor_r_tc_reject_rate = cor_rain_wl_con_wl_res$tc$reject_rate,
  cor_r_non_tc = cor_rain_wl_con_wl_res$non_tc$tau,
  cor_r_non_tc_reject_rate = cor_rain_wl_con_wl_res$non_tc$reject_rate)


cor_24hr_rain_ntr_con_24hr_rain <- list()
cor_24hr_rain_ntr_con_ntr <- list()
cor_24hr_rain_wl_con_24hr_rain <- list()
cor_24hr_rain_wl_con_wl <- list()

totals_fnc = function(Data, Window_Width_Sum){
  cs <- cumsum(Data)
  n  <- length(Data)
  
  Sum <- c(rep(NA, Window_Width_Sum - 1),
           cs[Window_Width_Sum:n] - c(0, cs[1:(n - Window_Width_Sum)]))
  
  return(list(Totals = Sum))
}


# Kendall correlation bootstrap con 24hr rainfall with ntr  ─────────────────────────────────────────────
library(energy)

# Pre-compute quantities used repeatedly
sim_rain_24hr_sums <- apply(sim_rain, 1, function(row) {
  max(totals_fnc(row, Window_Width_Sum = 24)$Totals, na.rm = TRUE)
})

sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
is_tc          <- aorc_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")


cor_24hr_rain_ntr_con_24hr_rain = list()
for (type in c("tc", "non_tc")) {
  
  obs_df      <- df_samples[df_samples$type_con_rain==type,c('rain_con_rain','ntr_con_rain')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_24hr_sums, wl$sim_ntr_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_24hr_sums, sim_ntr_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_24hr_rain_ntr_con_24hr_rain[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
}


res_cor_24hr_rain_ntr_con_24hr_rain_list[[i]] <- data.frame(
  cor_r_tc = cor_24hr_rain_ntr_con_24hr_rain$tc$tau,
  cor_r_tc_reject_rate = cor_24hr_rain_ntr_con_24hr_rain$tc$reject_rate,
  cor_r_non_tc = cor_24hr_rain_ntr_con_24hr_rain$non_tc$tau,
  cor_r_non_tc_reject_rate = cor_24hr_rain_ntr_con_24hr_rain$non_tc$reject_rate)


# Kendall correlation bootstrap 24hr rainfall con ntr 

sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")


cor_24hr_rain_ntr_con_ntr = list()
for (type in c("tc", "non_tc")) {
  obs_df      <- df_samples[df_samples$type_con_ntr==type,c('rain_con_ntr','ntr_con_ntr')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_24hr_sums, wl$sim_ntr_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_24hr_sums, sim_ntr_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_24hr_rain_ntr_con_ntr[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
  
}

res_cor_24hr_rain_ntr_con_ntr_list[[i]] <- data.frame(
  cor_ntr_tc = cor_24hr_rain_ntr_con_ntr$tc$tau,
  cor_ntr_tc_reject_rate = cor_24hr_rain_ntr_con_ntr$tc$reject_rate,
  cor_ntr_non_tc = cor_24hr_rain_ntr_con_ntr$non_tc$tau,
  cor_ntr_non_tc_reject_rate = cor_24hr_rain_ntr_con_ntr$non_tc$reject_rate)



# Kendall correlation bootstrap con 24hr rainfall with wl  ─────────────────────────────────────────────
library(energy)

# Pre-compute quantities used repeatedly
sim_rain_24hr_sums <- apply(sim_rain, 1, function(row) {
  max(totals_fnc(row, Window_Width_Sum = 24)$Totals, na.rm = TRUE)
})

sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
is_tc          <- aorc_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")


cor_24hr_rain_wl_con_24hr_rain = list()
for (type in c("tc", "non_tc")) {
  
  obs_df      <- df_samples[df_samples$type_con_rain==type,c('rain_con_rain','ntr_con_rain')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_24hr_sums, wl$sim_wl_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_24hr_sums, sim_wl_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_24hr_rain_wl_con_24hr_rain[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
}

res_cor_24hr_rain_wl_con_24hr_rain_list[[i]] <- data.frame(
  cor_r_tc = cor_24hr_rain_wl_con_24hr_rain$tc$tau,
  cor_r_tc_reject_rate = cor_24hr_rain_wl_con_24hr_rain$tc$reject_rate,
  cor_r_non_tc = cor_24hr_rain_wl_con_24hr_rain$non_tc$tau,
  cor_r_non_tc_reject_rate = cor_24hr_rain_wl_con_24hr_rain$non_tc$reject_rate)

# Kendall correlation bootstrap 24hr rainfall con wl

sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")


cor_24hr_rain_wl_con_wl = list()
for (type in c("tc", "non_tc")) {
  obs_df      <- df_samples[df_samples$type_con_ntr==type,c('rain_con_ntr','ntr_con_ntr')]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  
  # Simulated
  bt <- bootstrap_kendall(sim_rows, sim_rain_24hr_sums, wl$sim_wl_trim,
                          n = nrow(obs_df))
  
  #Energy Distance / E-statistic
  obs_mat_num <- as.matrix(obs_df)
  sim_df = data.frame(sim_rain_24hr_sums, sim_wl_max)
  sim_df = sim_df[sim_rows,]
  sim_mat_num <- as.matrix(sim_df)
  
  edist_test <- eqdist.etest(rbind(obs_mat_num, sim_mat_num),
                             sizes = c(nrow(obs_mat_num), nrow(sim_mat_num)),
                             R = 999)  # bootstrap replicates
  
  cor_24hr_rain_wl_con_wl[[type]] <- list(
    tau            = bt$tau,
    reject_rate = bt$reject_rate,
    edist_test_val = edist_test$estimate,
    edist_test_p_val = edist_test$p.value)
  
  
}

res_cor_24hr_rain_wl_con_wl_list[[i]] <- data.frame(
  cor_wl_tc = cor_24hr_rain_wl_con_wl$tc$tau,
  cor_wl_tc_reject_rate = cor_24hr_rain_wl_con_wl$tc$reject_rate,
  cor_wl_non_tc = cor_24hr_rain_wl_con_wl$non_tc$tau,
  cor_wl_non_tc_reject_rate = cor_24hr_rain_wl_con_wl$non_tc$reject_rate)

##Depth Duration curves
totals_fnc = function(Data, Window_Width_Sum){
  cs <- cumsum(Data)
  n  <- length(Data)
  
  Sum <- c(rep(NA, Window_Width_Sum - 1),
           cs[Window_Width_Sum:n] - c(0, cs[1:(n - Window_Width_Sum)]))
  
  return(list(Totals = Sum))
}

event_total_per_duration <- function(precip, start, end, dur){
  totals = array(0,dim=c(length(start),48))
  for(i in 1:length(start)){
    for(j in 1:min(dur[i],48)){
      totals[i,j] = max(totals_fnc(precip[(start[i]):(end[i])],Window_Width_Sum=j)$Totals,na.rm=T)
    }
    if(dur[i]<48){ 
      totals[i,((dur[i]+1):48)] = totals[i,dur[i]]
    }
  }
  return(list(Totals = totals))
}


aorc_rainfall_events_tc <- aorc_events[match(aorc_event_data$ids[which(aorc_event_data$types == 'tc')],aorc_events$event_24_id), ]

aorc_rainfall_events_non_tc <- aorc_events[match(aorc_event_data$ids[which(aorc_event_data$types == 'non_tc')],aorc_events$event_24_id), ]

tots_tc = event_total_per_duration(precip = aorc_precip_df$precip, 
                                   start = aorc_rainfall_events_tc$start,
                                   end = aorc_rainfall_events_tc$end, 
                                   dur=aorc_rainfall_events_tc$d)$Totals



tots_non_tc = event_total_per_duration(precip = aorc_precip_df$precip, 
                                       start = aorc_rainfall_events_non_tc$start, 
                                       end = aorc_rainfall_events_non_tc$end, 
                                       dur=aorc_rainfall_events_non_tc$d)$Totals

sim_tots = event_total_per_duration(precip = as.numeric(t(sim_rain)), 
                                    start = seq(1,145*10000,145), 
                                    end = seq(145,145*10000,145), 
                                    dur= rep(145,nrow(sim_rain)))$Totals


rps_bootstrap <- function(n, sim_rows, sim_rain_sums, perc = 0.5, n_boot=100){
  
  # Bootstrap storage
  gpd_rp_1  <- numeric(n_boot)
  gpd_rp_2  <- numeric(n_boot)
  gpd_rp_5  <- numeric(n_boot)
  gpd_rp_10  <- numeric(n_boot)
  gpd_rp_25  <- numeric(n_boot)
  gpd_rp_50  <- numeric(n_boot)
  gpd_rp_100  <- numeric(n_boot)
  gpd_rp_200  <- numeric(n_boot)
  gpd_rp_500  <- numeric(n_boot)
  gpd_rp_1000  <- numeric(n_boot)
  sim_sorted <- vector("list", n_boot)
  
  for(i in seq_len(n_boot)){
    
    samp            <- sample(length(sim_rows), n, replace=TRUE)
    s               <- sort(sim_rain_sums[sim_rows][samp])  # bootstrap sample
    sim_sorted[[i]] <- s
    #print(summary(s))
    gpd      <- GPD_Fit(Data=s, Data_Full=s, Thres=min(s), Method="Solari")  # fit on s
    gpd$Rate <- n / 44
    
    gpd_rp_1[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((1  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_2[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((2  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_5[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((5  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_10[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_25[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((25  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_50[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_100[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_200[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((200  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_500[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((500  * gpd$Rate)^gpd$xi - 1)
    gpd_rp_1000[i] <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((1000  * gpd$Rate)^gpd$xi - 1)
  }
  
  list(
    "rp1_50"  = quantile(gpd_rp_1, perc),
    "rp2_50"  = quantile(gpd_rp_2, perc),
    "rp5_50"  = quantile(gpd_rp_5, perc),
    "rp10_50"  = quantile(gpd_rp_10, perc),
    "rp25_50"  = quantile(gpd_rp_25, perc),
    "rp50_50"  = quantile(gpd_rp_50, perc),
    "rp100_50"  = quantile(gpd_rp_100, perc),
    "rp200_50"  = quantile(gpd_rp_200, perc),
    "rp500_50"  = quantile(gpd_rp_500, perc),
    "rp1000_50"  = quantile(gpd_rp_1000, perc)
  )
}

# Pre-compute quantities used repeatedly
sim_ddf_sums  <- sim_tots
is_tc          <- aorc_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")

# ── Rain marginal bootstrap ────────────────────────────────────────────────────
ddf_res <- list()
for (type in c("non_tc", "tc")) {
  obs      <- if (type =="tc"){tots_tc} else{tots_non_tc}
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- nrow(obs)
  rp       <- rp_axis(n)
  
  for (dur in 1:48) {
    
    # Observed GPD (on actual data, outside bootstrap loop)
    gpd_obs      <- GPD_Fit(Data=obs[,dur], Data_Full=obs[,dur], Thres=min(obs[,dur]), Method="Solari")
    gpd_obs$Rate <- n / 44
    obs_gpd_rp   <- gpd_obs$Threshold + (gpd_obs$sigma / gpd_obs$xi) * 
      ((c(1,2,5,10,25,50,100,200,500,1000) * gpd_obs$Rate)^gpd_obs$xi - 1)
    #print(summary(sim_ddf_sums[,dur]))
    
    bt <- rps_bootstrap(n=n, sim_rows, sim_ddf_sums[, dur])
    
    ddf_res[[type]][[paste0("dur_",dur)]] <- list(
      
      obs_rp1   = obs_gpd_rp[1],
      obs_rp2   = obs_gpd_rp[2],
      obs_rp5   = obs_gpd_rp[3], 
      obs_rp10  = obs_gpd_rp[4],
      obs_rp25  = obs_gpd_rp[5],
      obs_rp50  = obs_gpd_rp[6],
      obs_rp100 = obs_gpd_rp[7],  
      obs_rp200 = obs_gpd_rp[8],  
      obs_rp500 = obs_gpd_rp[9],  
      obs_rp1000 = obs_gpd_rp[10],  
      
      sim_rp1  = bt$rp1_50,
      sim_rp2  = bt$rp2_50,
      sim_rp5  = bt$rp5_50,
      sim_rp10  = bt$rp10_50,
      sim_rp25  = bt$rp25_50,
      sim_rp50  = bt$rp50_50,
      sim_rp100 = bt$rp100_50,
      sim_rp200 = bt$rp200_50,
      sim_rp500 = bt$rp500_50,
      sim_rp1000 = bt$rp1000_50
    )
  }
}

# inside your iteration loop:
res_ddf_list_tc[[i]] <- do.call(rbind, lapply(names(ddf_res$tc), function(d) {
  data.frame(aorc=i, dur=as.integer(gsub("dur_","",d)),
             as.data.frame(ddf_res$tc[[d]]))
}))

res_ddf_list_non_tc[[i]] <- do.call(rbind, lapply(names(ddf_res$non_tc), function(d) {
  data.frame(aorc=i, dur=as.integer(gsub("dur_","",d)),
             as.data.frame(ddf_res$non_tc[[d]]))
}))

sim_ddf_sums  <- sim_tots
rows_tc        <- which(df_sim$sample == "r_tc")
rows_non_tc    <- which(df_sim$sample == "r_non_tc")
rows_r_all     <- c(rows_tc, rows_non_tc)

n_tc     <- nrow(tots_tc)
n_non_tc <- nrow(tots_non_tc)
n_all    <- n_tc + n_non_tc

# Pool observed
tots_all <- rbind(tots_tc, tots_non_tc)

ddf_res <- list()

for (dur in 1:48) {
  
  # Observed — pooled
  gpd_obs <- GPD_Fit(Data=tots_all[,dur], Data_Full=tots_all[,dur], 
                     Thres=min(tots_all[,dur]), Method="Solari")
  gpd_obs$Rate <- n_all / 44
  obs_gpd_rp <- gpd_obs$Threshold + (gpd_obs$sigma / gpd_obs$xi) * 
    ((c(1,2,5,10,25,50,100,200,500,1000) * gpd_obs$Rate)^gpd_obs$xi - 1)
  
  # Simulated — pooled r samples, rate = n_all / 44
  bt <- rps_bootstrap(n=n_all, sim_rows=rows_r_all, sim_rain_sums=sim_ddf_sums[, dur])
  
  ddf_res[[paste0("dur_",dur)]] <- list(
    obs_rp1    = obs_gpd_rp[1],
    obs_rp2    = obs_gpd_rp[2],
    obs_rp5    = obs_gpd_rp[3],
    obs_rp10   = obs_gpd_rp[4],
    obs_rp25   = obs_gpd_rp[5],
    obs_rp50   = obs_gpd_rp[6],
    obs_rp100  = obs_gpd_rp[7],
    obs_rp200  = obs_gpd_rp[8],
    obs_rp500  = obs_gpd_rp[9],
    obs_rp1000 = obs_gpd_rp[10],
    
    sim_rp1    = bt$rp1_50,
    sim_rp2    = bt$rp2_50,
    sim_rp5    = bt$rp5_50,
    sim_rp10   = bt$rp10_50,
    sim_rp25   = bt$rp25_50,
    sim_rp50   = bt$rp50_50,
    sim_rp100  = bt$rp100_50,
    sim_rp200  = bt$rp200_50,
    sim_rp500  = bt$rp500_50,
    sim_rp1000 = bt$rp1000_50
  )
}

# Store results — single pooled DDF per AORC point
res_ddf_all_con_r_list[[i]] <- do.call(rbind, lapply(names(ddf_res), function(d) {
  data.frame(aorc=i, dur=as.integer(gsub("dur_","",d)),
             as.data.frame(ddf_res[[d]]))
}))

}

#cor_con_ntr_tc_edist_test_val = cor_con_ntr_res$tc$edist_test_val,
#cor_con_ntr_tc_edist_test_p_val = cor_con_ntr_res$tc$edist_test_p_val,
#cor_con_ntr_non_tc_edist_test_val = cor_con_ntr_res$non_tc$edist_test_val,
#cor_con_ntr_non_tc_edist_test_p_val = cor_con_ntr_res$non_tc$edist_test_p_val

res_cor_con_rain <- do.call(rbind, res_cor_con_rain_list)
res_cor_con_ntr  <- do.call(rbind, res_cor_con_ntr_list)
res_cor_rain_wl_con_rain  <- do.call(rbind,res_cor_rain_wl_con_rain_list)
res_cor_rain_wl_con_wl    <- do.call(rbind,res_cor_rain_wl_con_wl_list)

res_cor_24hr_rain_ntr_con_ntr        <- do.call(rbind, res_cor_24hr_rain_ntr_con_ntr_list)
res_cor_24hr_rain_ntr_con_24hr_rain  <- do.call(rbind, res_cor_24hr_rain_ntr_con_24hr_rain_list)
res_cor_24hr_rain_wl_con_wl          <- do.call(rbind, res_cor_24hr_rain_wl_con_wl_list)
res_cor_24hr_rain_wl_con_24hr_rain   <- do.call(rbind, res_cor_24hr_rain_wl_con_24hr_rain_list)

res_rain <- do.call(rbind, res_rain_list)
res_rain_df_tc <- do.call(rbind, res_rain_df_tc_list)
res_rain_df_non_tc <- do.call(rbind, res_rain_df_non_tc_list)

res_vol <- do.call(rbind, res_vol_list)
res_vol_df_tc <- do.call(rbind, res_vol_df_tc_list)
res_vol_df_non_tc <- do.call(rbind, res_vol_df_non_tc_list)


# At the end, combine across iterations:
res_ddf_tc     <- do.call(rbind, res_ddf_list_tc)
res_ddf_non_tc <- do.call(rbind, res_ddf_list_non_tc)
res_ddf_all_con_r <- do.call(rbind, res_ddf_all_con_r_list)

write.csv(res_cor_con_rain, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_con_rain_aorc_subset.csv',sep=""))
write.csv(res_cor_con_ntr, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_con_ntr_aorc_subset.csv',sep=""))
write.csv(res_cor_rain_wl_con_rain, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_rain_wl_con_rain_aorc_subset.csv',sep=""))
write.csv(res_cor_rain_wl_con_wl, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_rain_wl_con_wl_aorc_subset.csv',sep=""))

write.csv(res_cor_24hr_rain_ntr_con_ntr, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_ntr_aorc_subset.csv',sep=""))
write.csv(res_cor_24hr_rain_ntr_con_24hr_rain, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_24hr_rain_aorc_subset.csv',sep=""))
write.csv(res_cor_24hr_rain_wl_con_wl, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_wl_con_wl_aorc_subset.csv',sep=""))
write.csv(res_cor_24hr_rain_wl_con_24hr_rain , paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_wl_con_24hr_rain_aorc_subset.csv',sep=""))


write.csv(res_vol, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_aorc_subset.csv',sep=""))
write.csv(res_vol_df_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_df_tc_aorc_subset.csv',sep=""))
write.csv(res_vol_df_non_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_df_non_tc_aorc_subset.csv',sep=""))

write.csv(res_rain, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_aorc_subset.csv',sep=""))
write.csv(res_rain_df_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_tc_aorc_subset.csv',sep=""))
write.csv(res_rain_df_non_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_non_tc_aorc_subset.csv',sep=""))

write.csv(res_ddf_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ddf_tc_aorc_subset.csv',sep=""))
write.csv(res_ddf_non_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ddf_non_tc_aorc_subset.csv',sep=""))

write.csv(res_ddf_all_con_r, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ddf_all_con_r_aorc_subset.csv',sep=""))


res_ntr_list = list()
res_wl_list = list()
res_ntr_df_tc_list = list()
res_wl_df_tc_list = list()
res_ntr_df_non_tc_list = list()
res_wl_df_non_tc_list = list()

# ── ntr marginal bootstrap ────────────────────────────────────────────────────
# Pre-compute quantities used repeatedly
sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
is_tc          <- ntr_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")

ntr_res <- list()
for (type in c("non_tc", "tc")) {
  
  obs      <- cora_ntr_df$ntr[ntr_event_data$ids][if (type == "tc") is_tc else !is_tc]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- length(obs)
  rp       <- rp_axis(n)
  
  s               <- sort(obs)
  obs_rp10        <- approx(rp, s, xout = 10)$y
  obs_rp40         <- approx(rp, s, xout = 40)$y
  print(summary(s))
  print("Here")
  gpd <- GPD_Fit(Data=s, Data_Full=cora_ntr_df$ntr, Thres=min(s), Method="Solari")
  gpd$Rate <- n/44
  print("Here_1")
  #Extract ten year event
  obs_gpd_rp10  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp40  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((40*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp50  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp100 <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100*gpd$Rate)^gpd$xi - 1)  
  
  obs_gpd_threshold  <- gpd$Threshold 
  obs_gpd_xi  <- gpd$xi
  obs_gpd_sigma <- gpd$sigma
  obs_ar2 <- AR2(c(obs_gpd_xi,obs_gpd_sigma,obs_gpd_threshold),obs)
  
  obs_perc = ((1:n)-0.5)/n
  
  bt <- bootstrap_rain(obs, sim_rows, sim_rain_sums=sim_ntr_max, data=cora_ntr_df$ntr, perc=obs_perc)
  
  ntr_res[[type]] <- list(
    obs = obs,
    obs_rp10           = obs_rp10,
    obs_rp40           = obs_rp40,
    obs_gpd_rp10       = obs_gpd_rp10, 
    obs_gpd_rp40       = obs_gpd_rp40, 
    obs_gpd_rp50       = obs_gpd_rp50,  
    obs_gpd_rp100      = obs_gpd_rp100,
    obs_gpd_threshold  = obs_gpd_threshold, 
    obs_gpd_xi         = obs_gpd_xi,
    obs_gpd_sigma      = obs_gpd_sigma, 
    obs_ar2            = obs_ar2, 
    
    rp10_025      = quantile(bt$rp10,0.025),
    rp10_05      = quantile(bt$rp10,0.05),
    rp10_25      = quantile(bt$rp10,0.25),
    rp10_50      = quantile(bt$rp10,0.5),
    rp10_75      = quantile(bt$rp10,0.75),
    rp10_95      = quantile(bt$rp10,0.95),
    rp10_975      = quantile(bt$rp10,0.975),
    rp40_025      = quantile(bt$rp40,0.025),
    rp40_05      = quantile(bt$rp40,0.05),
    rp40_25      = quantile(bt$rp40,0.25),
    rp40_50      = quantile(bt$rp40,0.5),
    rp40_75      = quantile(bt$rp40,0.75),
    rp40_95      = quantile(bt$rp40,0.95),
    rp40_975      = quantile(bt$rp40,0.975),
    gpd_rp10_025       = quantile(bt$gpd_rp10,0.025), 
    gpd_rp10_05       = quantile(bt$gpd_rp10,0.05), 
    gpd_rp10_25       = quantile(bt$gpd_rp10,0.25), 
    gpd_rp10_50       = quantile(bt$gpd_rp10,0.5),
    gpd_rp10_75       = quantile(bt$gpd_rp10,0.75),
    gpd_rp10_95       = quantile(bt$gpd_rp10,0.95), 
    gpd_rp10_975       = quantile(bt$gpd_rp10,0.975), 
    gpd_rp40_025       = quantile(bt$gpd_rp40,0.025),
    gpd_rp40_05       = quantile(bt$gpd_rp40,0.05),
    gpd_rp40_25       = quantile(bt$gpd_rp40,0.25),
    gpd_rp40_50       = quantile(bt$gpd_rp40,0.50),
    gpd_rp40_75       = quantile(bt$gpd_rp40,0.75),
    gpd_rp40_95       = quantile(bt$gpd_rp40,0.95),
    gpd_rp40_975       = quantile(bt$gpd_rp40,0.975),
    gpd_rp50_025       = quantile(bt$gpd_rp50,0.025),
    gpd_rp50_05       = quantile(bt$gpd_rp50,0.05), 
    gpd_rp50_25       = quantile(bt$gpd_rp50,0.25), 
    gpd_rp50_50       = quantile(bt$gpd_rp50,0.5),
    gpd_rp50_75       = quantile(bt$gpd_rp50,0.75),
    gpd_rp50_95       = quantile(bt$gpd_rp50,0.95),
    gpd_rp50_975       = quantile(bt$gpd_rp50,0.975),
    gpd_rp100_025      = quantile(bt$gpd_rp100,0.025),
    gpd_rp100_05      = quantile(bt$gpd_rp100,0.05),
    gpd_rp100_25      = quantile(bt$gpd_rp100,0.25),
    gpd_rp100_50      = quantile(bt$gpd_rp100,0.50),
    gpd_rp100_75      = quantile(bt$gpd_rp100,0.75),
    gpd_rp100_95      = quantile(bt$gpd_rp100,0.95),
    gpd_rp100_975      = quantile(bt$gpd_rp100,0.975),
    gpd_threshold_025  = quantile(bt$gpd_threshold, 0.025),
    gpd_threshold_05  = quantile(bt$gpd_threshold, 0.05),
    gpd_threshold_25  = quantile(bt$gpd_threshold, 0.25),
    gpd_threshold_50  = quantile(bt$gpd_threshold, 0.50),
    gpd_threshold_75  = quantile(bt$gpd_threshold, 0.75),
    gpd_threshold_95  = quantile(bt$gpd_threshold, 0.95),
    gpd_threshold_975  = quantile(bt$gpd_threshold, 0.975),
    gpd_xi_025         = quantile(bt$gpd_xi,0.025),
    gpd_xi_05         = quantile(bt$gpd_xi,0.05),
    gpd_xi_25         = quantile(bt$gpd_xi,0.25),
    gpd_xi_50         = quantile(bt$gpd_xi,0.5),
    gpd_xi_75         = quantile(bt$gpd_xi,0.75),
    gpd_xi_95         = quantile(bt$gpd_xi,0.95),
    gpd_xi_975         = quantile(bt$gpd_xi,0.975),
    gpd_sigma_025      = quantile(bt$gpd_sigma, 0.025),
    gpd_sigma_05      = quantile(bt$gpd_sigma, 0.05),
    gpd_sigma_25      = quantile(bt$gpd_sigma, 0.25),
    gpd_sigma_50      = quantile(bt$gpd_sigma, 0.5),
    gpd_sigma_75      = quantile(bt$gpd_sigma, 0.75),
    gpd_sigma_95      = quantile(bt$gpd_sigma, 0.95),
    gpd_sigma_975      = quantile(bt$gpd_sigma, 0.975),
    ar2_025            = quantile(bt$ar2, 0.025),
    ar2_05            = quantile(bt$ar2, 0.05),
    ar2_25            = quantile(bt$ar2, 0.25),
    ar2_50            = quantile(bt$ar2, 0.5),
    ar2_75            = quantile(bt$ar2, 0.75),
    ar2_95            = quantile(bt$ar2, 0.95),
    ar2_975            = quantile(bt$ar2, 0.975),
    p_025           =  bt$p_025,
    p_05            =  bt$p_05,
    p_50            =  bt$p_50,
    p_95            =  bt$p_95,
    p_975           =  bt$p_975,
    p_val_reject = bt$p_val_reject
  )
}

res_ntr_list <- data.frame(
  
  obs_rp10_tc           =  ntr_res$tc$obs_rp10,
  obs_rp40_tc           =  ntr_res$tc$obs_rp40,
  obs_gpd_rp10_tc       =  ntr_res$tc$obs_gpd_rp10, 
  obs_gpd_rp40_tc       =  ntr_res$tc$obs_gpd_rp40, 
  obs_gpd_rp50_tc       =  ntr_res$tc$obs_gpd_rp50,  
  obs_gpd_rp100_tc      =  ntr_res$tc$obs_gpd_rp100,
  obs_gpd_threshold_tc  =  ntr_res$tc$obs_gpd_threshold, 
  obs_gpd_xi_tc         =  ntr_res$tc$obs_gpd_xi,
  obs_gpd_sigma_tc      =  ntr_res$tc$obs_gpd_sigma, 
  obs_ar2_tc            =  ntr_res$tc$obs_ar2, 
 
  rp10_025_tc  =  ntr_res$tc$rp10_025, 
  rp10_05_tc  =  ntr_res$tc$rp10_05,          
  rp10_25_tc  =  ntr_res$tc$rp10_25,           
  rp10_50_tc =   ntr_res$tc$rp10_50,          
  rp10_75_tc  =  ntr_res$tc$rp10_75,      
  rp10_95_tc  =  ntr_res$tc$rp10_95,
  rp10_975_tc  =  ntr_res$tc$rp10_975,
  
  rp40_025_tc  =  ntr_res$tc$rp40_025, 
  rp40_05_tc  =  ntr_res$tc$rp40_05, 
  rp40_25_tc  =  ntr_res$tc$rp40_25,         
  rp40_50_tc  =  ntr_res$tc$rp40_50,      
  rp40_75_tc  =  ntr_res$tc$rp40_75,      
  rp40_95_tc  =  ntr_res$tc$rp40_95, 
  rp40_975_tc  =  ntr_res$tc$rp40_975, 
  
  gpd_rp10_025_tc  =  ntr_res$tc$gpd_rp10_025, 
  gpd_rp10_05_tc  =  ntr_res$tc$gpd_rp10_05,        
  gpd_rp10_25_tc  =  ntr_res$tc$gpd_rp10_25,        
  gpd_rp10_50_tc  =  ntr_res$tc$gpd_rp10_50,      
  gpd_rp10_75_tc  =  ntr_res$tc$gpd_rp10_75,       
  gpd_rp10_95_tc  =  ntr_res$tc$gpd_rp10_95,        
  gpd_rp10_975_tc  =  ntr_res$tc$gpd_rp10_975, 
  
  gpd_rp40_025_tc  =  ntr_res$tc$gpd_rp40_025,       
  gpd_rp40_05_tc  =  ntr_res$tc$gpd_rp40_05,       
  gpd_rp40_25_tc =   ntr_res$tc$gpd_rp40_25,      
  gpd_rp40_50_tc   = ntr_res$tc$gpd_rp40_50,       
  gpd_rp40_75_tc   = ntr_res$tc$gpd_rp40_75,       
  gpd_rp40_95_tc   = ntr_res$tc$gpd_rp40_95, 
  gpd_rp40_975_tc   = ntr_res$tc$gpd_rp40_975, 
  
  gpd_rp50_025_tc  =  ntr_res$tc$gpd_rp50_025,   
  gpd_rp50_05_tc  =  ntr_res$tc$gpd_rp50_05,      
  gpd_rp50_25_tc  =  ntr_res$tc$gpd_rp50_25,       
  gpd_rp50_50_tc  =  ntr_res$tc$gpd_rp50_50,       
  gpd_rp50_75_tc  =  ntr_res$tc$gpd_rp50_75,       
  gpd_rp50_95_tc  =  ntr_res$tc$gpd_rp50_95, 
  gpd_rp50_975_tc  =  ntr_res$tc$gpd_rp50_975,   
  
  gpd_rp100_025_tc  =  ntr_res$tc$gpd_rp100_025, 
  gpd_rp100_05_tc  =  ntr_res$tc$gpd_rp100_05,      
  gpd_rp100_25_tc  =  ntr_res$tc$gpd_rp100_25,      
  gpd_rp100_50_tc  =  ntr_res$tc$gpd_rp100_50,      
  gpd_rp100_75_tc  =  ntr_res$tc$gpd_rp100_75,      
  gpd_rp100_95_tc  =  ntr_res$tc$gpd_rp100_95,
  gpd_rp100_975_tc  =  ntr_res$tc$gpd_rp100_975, 
  
  gpd_threshold_025_tc  =  ntr_res$tc$gpd_threshold_025,  
  gpd_threshold_05_tc  =  ntr_res$tc$gpd_threshold_05,  
  gpd_threshold_25_tc  =  ntr_res$tc$gpd_threshold_25, 
  gpd_threshold_50_tc  =  ntr_res$tc$gpd_threshold_50,  
  gpd_threshold_75_tc  =  ntr_res$tc$gpd_threshold_75,  
  gpd_threshold_95_tc  =  ntr_res$tc$gpd_threshold_95,  
  gpd_threshold_975_tc  =  ntr_res$tc$gpd_threshold_975,  
  
  gpd_xi_025_tc  =  ntr_res$tc$gpd_xi_025,     
  gpd_xi_05_tc  =  ntr_res$tc$gpd_xi_05,         
  gpd_xi_25_tc  =  ntr_res$tc$gpd_xi_25,         
  gpd_xi_50_tc  =  ntr_res$tc$gpd_xi_50,         
  gpd_xi_75_tc  =  ntr_res$tc$gpd_xi_75,        
  gpd_xi_95_tc  =  ntr_res$tc$gpd_xi_95,     
  gpd_xi_975_tc  =  ntr_res$tc$gpd_xi_975,     
  
  gpd_sigma_025_tc  =  ntr_res$tc$gpd_sigma_025,  
  gpd_sigma_05_tc  =  ntr_res$tc$gpd_sigma_05,     
  gpd_sigma_25_tc  =  ntr_res$tc$gpd_sigma_25,     
  gpd_sigma_50_tc  =  ntr_res$tc$gpd_sigma_50,     
  gpd_sigma_75_tc  =  ntr_res$tc$gpd_sigma_75,     
  gpd_sigma_95_tc  =  ntr_res$tc$gpd_sigma_95,  
  gpd_sigma_975_tc  =  ntr_res$tc$gpd_sigma_975,  
  
  ar2_025_tc    =  ntr_res$tc$ar2_025, 
  ar2_05_tc    =  ntr_res$tc$ar2_05,            
  ar2_25_tc    =  ntr_res$tc$ar2_25,            
  ar2_50_tc    =  ntr_res$tc$ar2_50,            
  ar2_75_tc    =  ntr_res$tc$ar2_75,            
  ar2_95_tc    =  ntr_res$tc$ar2_95, 
  ar2_975_tc    =  ntr_res$tc$ar2_975, 
  
  p_val_reject_tc  =  ntr_res$tc$p_val_reject,
  
  obs_rp10_non_tc          =  ntr_res$non_tc$obs_rp10,
  obs_rp40_non_tc          =  ntr_res$non_tc$obs_rp40,
  obs_gpd_rp10_non_tc      =  ntr_res$non_tc$obs_gpd_rp10, 
  obs_gpd_rp40_non_tc      =  ntr_res$non_tc$obs_gpd_rp40, 
  obs_gpd_rp50_non_tc      =  ntr_res$non_tc$obs_gpd_rp50,  
  obs_gpd_rp100_non_tc     =  ntr_res$non_tc$obs_gpd_rp100,
  obs_gpd_threshold_non_tc  =  ntr_res$non_tc$obs_gpd_threshold, 
  obs_gpd_xi_non_tc         =  ntr_res$non_tc$obs_gpd_xi,
  obs_gpd_sigma_non_tc      =  ntr_res$non_tc$obs_gpd_sigma, 
  obs_ar2_non_tc            =  ntr_res$non_tc$obs_ar2, 
  
  rp10_025_non_tc  =  ntr_res$non_tc$rp10_025,     
  rp10_05_non_tc  =  ntr_res$non_tc$rp10_05,          
  rp10_25_non_tc  =  ntr_res$non_tc$rp10_25,           
  rp10_50_non_tc  =  ntr_res$non_tc$rp10_50,          
  rp10_75_non_tc  =  ntr_res$non_tc$rp10_75,      
  rp10_95_non_tc  =  ntr_res$non_tc$rp10_95, 
  rp10_975_non_tc  =  ntr_res$non_tc$rp10_975,     
  
  rp40_025_non_tc  =   ntr_res$non_tc$rp40_025, 
  rp40_05_non_tc  =   ntr_res$non_tc$rp40_05, 
  rp40_25_non_tc  =   ntr_res$non_tc$rp40_25,         
  rp40_50_non_tc  =   ntr_res$non_tc$rp40_50,      
  rp40_75_non_tc  =   ntr_res$non_tc$rp40_75,      
  rp40_95_non_tc  =   ntr_res$non_tc$rp40_95, 
  rp40_975_non_tc  =   ntr_res$non_tc$rp40_975, 
  
  gpd_rp10_025_non_tc  =  ntr_res$non_tc$gpd_rp10_025,
  gpd_rp10_05_non_tc  =  ntr_res$non_tc$gpd_rp10_05,
  gpd_rp10_25_non_tc  =  ntr_res$non_tc$gpd_rp10_25,        
  gpd_rp10_50_non_tc  =  ntr_res$non_tc$gpd_rp10_50,      
  gpd_rp10_75_non_tc  =  ntr_res$non_tc$gpd_rp10_75,       
  gpd_rp10_95_non_tc  =  ntr_res$non_tc$gpd_rp10_95,
  gpd_rp10_975_non_tc  =  ntr_res$non_tc$gpd_rp10_975,
  
  gpd_rp40_025_non_tc  =   ntr_res$non_tc$gpd_rp40_025,  
  gpd_rp40_05_non_tc  =   ntr_res$non_tc$gpd_rp40_05,       
  gpd_rp40_25_non_tc  =   ntr_res$non_tc$gpd_rp40_25,      
  gpd_rp40_50_non_tc  =   ntr_res$non_tc$gpd_rp40_50,       
  gpd_rp40_75_non_tc  =   ntr_res$non_tc$gpd_rp40_75,       
  gpd_rp40_95_non_tc  =   ntr_res$non_tc$gpd_rp40_95,  
  gpd_rp40_975_non_tc  =   ntr_res$non_tc$gpd_rp40_975,  
  
  gpd_rp50_025_non_tc  =   ntr_res$non_tc$gpd_rp50_025,
  gpd_rp50_05_non_tc  =   ntr_res$non_tc$gpd_rp50_05,      
  gpd_rp50_25_non_tc  =   ntr_res$non_tc$gpd_rp50_25,       
  gpd_rp50_50_non_tc  =   ntr_res$non_tc$gpd_rp50_50,       
  gpd_rp50_75_non_tc  =   ntr_res$non_tc$gpd_rp50_75,       
  gpd_rp50_95_non_tc  =   ntr_res$non_tc$gpd_rp50_95,
  gpd_rp50_975_non_tc  =   ntr_res$non_tc$gpd_rp50_975,
  
  gpd_rp100_025_non_tc  =   ntr_res$non_tc$gpd_rp100_025,     
  gpd_rp100_05_non_tc  =   ntr_res$non_tc$gpd_rp100_05,      
  gpd_rp100_25_non_tc  =   ntr_res$non_tc$gpd_rp100_25,      
  gpd_rp100_50_non_tc  =   ntr_res$non_tc$gpd_rp100_50,      
  gpd_rp100_75_non_tc  =   ntr_res$non_tc$gpd_rp100_75,      
  gpd_rp100_95_non_tc  =   ntr_res$non_tc$gpd_rp100_95,     
  gpd_rp100_975_non_tc  =   ntr_res$non_tc$gpd_rp100_975,     
  
  gpd_threshold_025_non_tc  =   ntr_res$non_tc$gpd_threshold_025,  
  gpd_threshold_05_non_tc  =   ntr_res$non_tc$gpd_threshold_05,  
  gpd_threshold_25_non_tc  =   ntr_res$non_tc$gpd_threshold_25, 
  gpd_threshold_50_non_tc  =   ntr_res$non_tc$gpd_threshold_50,  
  gpd_threshold_75_non_tc  =   ntr_res$non_tc$gpd_threshold_75,  
  gpd_threshold_95_non_tc  =   ntr_res$non_tc$gpd_threshold_95,  
  gpd_threshold_975_non_tc  =   ntr_res$non_tc$gpd_threshold_975,  
  
  gpd_xi_025_non_tc  =  ntr_res$non_tc$gpd_xi_025, 
  gpd_xi_05_non_tc  =  ntr_res$non_tc$gpd_xi_05,         
  gpd_xi_25_non_tc  =  ntr_res$non_tc$gpd_xi_25,         
  gpd_xi_50_non_tc  =  ntr_res$non_tc$gpd_xi_50,         
  gpd_xi_75_non_tc  =  ntr_res$non_tc$gpd_xi_75,        
  gpd_xi_95_non_tc  =  ntr_res$non_tc$gpd_xi_95,
  gpd_xi_975_non_tc  =  ntr_res$non_tc$gpd_xi_975, 
  
  gpd_sigma_025_non_tc  =   ntr_res$non_tc$gpd_sigma_025, 
  gpd_sigma_05_non_tc  =   ntr_res$non_tc$gpd_sigma_05,     
  gpd_sigma_25_non_tc  =   ntr_res$non_tc$gpd_sigma_25,     
  gpd_sigma_50_non_tc  =   ntr_res$non_tc$gpd_sigma_50,     
  gpd_sigma_75_non_tc  =   ntr_res$non_tc$gpd_sigma_75,     
  gpd_sigma_95_non_tc  =   ntr_res$non_tc$gpd_sigma_95, 
  gpd_sigma_975_non_tc  =   ntr_res$non_tc$gpd_sigma_975, 
  
  ar2_025_non_tc   =  ntr_res$non_tc$ar2_025,
  ar2_05_non_tc   =  ntr_res$non_tc$ar2_05,            
  ar2_25_non_tc   =  ntr_res$non_tc$ar2_25,            
  ar2_50_non_tc   =  ntr_res$non_tc$ar2_50,            
  ar2_75_non_tc   =  ntr_res$non_tc$ar2_75,            
  ar2_95_non_tc   =  ntr_res$non_tc$ar2_95,
  ar2_975_non_tc   =  ntr_res$non_tc$ar2_975,
  
  p_val_reject_non_tc = ntr_res$non_tc$p_val_reject)

res_ntr_df_tc_list = data.frame(obs_tc            =  ntr_res$tc$obs,
                                p025_tc           =  ntr_res$tc$p_025,
                                p05_tc            =  ntr_res$tc$p_05,
                                p50_tc            =  ntr_res$tc$p_50,
                                p95_tc            =  ntr_res$tc$p_95,
                                p975_tc           =  ntr_res$tc$p_975)

res_ntr_df_non_tc_list  = data.frame(obs_non_tc            =  ntr_res$non_tc$obs,
                                     p025_non_tc           =  ntr_res$non_tc$p_025,
                                     p05_non_tc            =  ntr_res$non_tc$p_05,
                                     p50_non_tc            =  ntr_res$non_tc$p_50,
                                     p95_non_tc            =  ntr_res$non_tc$p_95,
                                     p975_non_tc           =  ntr_res$non_tc$p_975)



sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
is_tc          <- ntr_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")

ntr_gpd_res <- list()
for (type in c("non_tc", "tc")) {
  
  obs      <- cora_ntr_df$ntr[ntr_event_data$ids][if (type == "tc") is_tc else !is_tc]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- length(obs)
  rp       <- rp_axis(n)
  
  s               <- sort(obs)
  gpd <- evm(s, th=min(s))
  gpd$rate <- n/44
  print("Gpd_obs")
  print(gpd)
  obs_rl_lower = numeric(100)
  obs_rl_est = numeric(100)
  obs_rl_upper = numeric(100)
  rp= 1:100
  for(i in 1:100){
    
    rl_out <- tryCatch({
      
      out <- rl(gpd,
                M = gpd$rate * rp[i],
                alpha = 0.05,
                ci.fit = TRUE)[[1]]
      
      c(est   = out[1],
        lower = out[2],
        upper = out[3])
      
    }, error = function(e) {
      
      if (grepl("M \\* rate must be > 1", e$message)) {
        c(est = NA, lower = NA, upper = NA)
      } else {
        stop(e)
      }
      
    })
    
    obs_rl_est[i]   <- rl_out["est"]
    obs_rl_lower[i] <- rl_out["lower"]
    obs_rl_upper[i] <- rl_out["upper"]
  }
  
  gpd <- evm(sim_ntr_max[sim_rows], th=min(sim_ntr_max[sim_rows])) 
  gpd$rate <- n/44
  print("Gpd_sim")
  print(gpd)
  sim_rl_lower = numeric(100)
  sim_rl_est = numeric(100)
  sim_rl_upper = numeric(100)
  rp= 1:100
  for(i in 1:100){
    
    rl_out <- tryCatch({
      
      out <- rl(gpd,
                M = gpd$rate * rp[i],
                alpha = 0.05,
                ci.fit = TRUE)[[1]]
      
      c(est   = out[1],
        lower = out[2],
        upper = out[3])
      
    }, error = function(e) {
      
      if (grepl("M \\* rate must be > 1", e$message)) {
        c(est = NA, lower = NA, upper = NA)
      } else {
        stop(e)
      }
      
    })
    
    sim_rl_est[i]   <- rl_out["est"]
    sim_rl_lower[i] <- rl_out["lower"]
    sim_rl_upper[i] <- rl_out["upper"]
  }
  
  ntr_gpd_res[[type]] <- list(
    obs_rl_lower = obs_rl_lower,
    obs_rl_est  = obs_rl_est,
    obs_rl_upper = obs_rl_upper,
    
    sim_rl_lower = sim_rl_lower,
    sim_rl_est = sim_rl_est,
    sim_rl_upper = sim_rl_upper
  )
}

res_ntr_gpd= data.frame(obs_rl_lower_tc    =  ntr_gpd_res$tc$obs_rl_lower,
                        obs_rl_est_tc      =  ntr_gpd_res$tc$obs_rl_est,
                        obs_rl_upper_tc    =  ntr_gpd_res$tc$obs_rl_upper,
                         
                        sim_rl_lower_tc       =  ntr_gpd_res$tc$sim_rl_lower,
                        sim_rl_est_tc         =  ntr_gpd_res$tc$sim_rl_est,
                        sim_rl_upper_tc       =  ntr_gpd_res$tc$sim_rl_upper,
                         
                        obs_rl_lower_non_tc     =  ntr_gpd_res$non_tc$obs_rl_lower,
                        obs_rl_est_non_tc       =  ntr_gpd_res$non_tc$obs_rl_est,
                        obs_rl_upper_non_tc     =  ntr_gpd_res$non_tc$obs_rl_upper,
                         
                        sim_rl_lower_non_tc     =  ntr_gpd_res$non_tc$sim_rl_lower,
                        sim_rl_est_non_tc       =  ntr_gpd_res$non_tc$sim_rl_est,
                        sim_rl_upper_non_tc     =  ntr_gpd_res$non_tc$sim_rl_upper)

write.csv(res_ntr_gpd, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_ntr_gpd.csv',sep=""))

#Obs gpd combined
obs = cora_ntr_df$ntr[ntr_event_data$ids]
n = length(obs)

#Probabilities for multinomial distribution
total_events = 5*44
p_ntr_tc = length(which(ntr_event_data$types == "tc"))/total_events
p_ntr_non_tc = length(which(ntr_event_data$types == "non_tc"))/total_events
probs = c(p_ntr_tc, p_ntr_non_tc)

obs_rl_lower = array(NA, dim=c(100,100))
obs_rl_est = array(NA, dim=c(100,100))
obs_rl_upper = array(NA, dim=c(100,100))

rp= 1:100

#Pull out indices of TC and non-tc conditioned 
idx_tc_obs  = which(ntr_event_data$types == "tc")
idx_non_tc_obs = which(ntr_event_data$types == "non_tc")

for(i in 1:100){ 
  #Sample from multinomial
  type = rmultinom(1, size=100*5, prob = probs)
  
  #Number of each type
  n_tc = type[1]
  n_non_tc = type[2]
  
  #Pull out indicies of sampled rows
  rows = c(sample(idx_tc_obs,n_tc,replace=T),
           sample(idx_non_tc_obs,n_non_tc,replace=T))
  
  #Fit GPD
  gpd <- evm(obs[rows], th=min(obs[rows])) 
  gpd$rate <- n/44
  
  #Estimate retunr levels
  for(j in 1:100){
    rl_out = rl(gpd, M = gpd$rate*rp[j], alpha = 0.05, ci.fit = TRUE)[[1]]
    obs_rl_lower[i,j] = rl_out[2]
    obs_rl_est[i,j] = rl_out[1]
    obs_rl_upper[i,j] = rl_out[3]
  }
} 

write.csv(obs_rl_lower, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_lower.csv',sep=""))
write.csv(obs_rl_est, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_est.csv',sep=""))
write.csv(obs_rl_upper, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_upper.csv',sep=""))



#Simulations

sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)

obs      <- cora_ntr_df$ntr[ntr_event_data$ids]
n        <- length(obs)

#Probabilities for multinomial distribution
total_events = 5*44
p_ntr_tc = length(which(ntr_event_data$types == "tc"))/total_events
p_ntr_non_tc = length(which(ntr_event_data$types == "non_tc"))/total_events
probs = c(p_ntr_tc, p_ntr_non_tc)
 
sim_rl_lower = array(NA, dim=c(100,100))
sim_rl_est = array(NA, dim=c(100,100))
sim_rl_upper = array(NA, dim=c(100,100))

rp= 1:100

#Pull out indicies of Tc and non-TC events
idx_tc_sim = which(df_sim$sample == "ntr_tc")
idx_non_tc_sim = which(df_sim$sample == "ntr_non_tc")     

for(i in 1:100){ 
 #Sample from multinomial
 type = rmultinom(1, size=100*5, prob = probs)
 
 #Pull out number of sampled Tc and non-TC events
 n_tc = type[1]
 n_non_tc = type[2]       
 
 #Get indicies of sampled rows
 rows = c(sample(idx_tc_sim,n_tc, replace=T),
          sample(idx_non_tc_sim,n_non_tc, replace=T))   
 
 #Fit GPD
 gpd <- evm(sim_ntr_max[rows], th=min(sim_ntr_max[rows])) 
 gpd$rate <- n/44
 
 #Estimate return levels
 for(j in 1:100){
    rl_out = rl(gpd, M = gpd$rate*rp[j], alpha = 0.05, ci.fit = TRUE)[[1]]
    sim_rl_lower[i,j] = rl_out[2]
    sim_rl_est[i,j] = rl_out[1]
    sim_rl_upper[i,j] = rl_out[3]
  }
} 

write.csv(sim_rl_lower, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_lower.csv',sep=""))
write.csv(sim_rl_est, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_est.csv',sep=""))
write.csv(sim_rl_upper, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_upper.csv',sep=""))


#Max observed wl within +/- 72 hours of ntr peak
obs = sapply(ntr_event_data$ids,function(x){max(cora_df$wl[(x-72):(x+72)])})
n = length(obs)

#Probabilities for multinomial distribution
total_events = 5*44
p_ntr_tc = length(which(ntr_event_data$types == "tc"))/total_events
p_ntr_non_tc = length(which(ntr_event_data$types == "non_tc"))/total_events
probs = c(p_ntr_tc, p_ntr_non_tc)

obs_rl_lower = array(NA, dim=c(100,100))
obs_rl_est = array(NA, dim=c(100,100))
obs_rl_upper = array(NA, dim=c(100,100))

rp= 1:100

#Pull out indices of TC and non-tc conditioned 
idx_tc_obs  = which(ntr_event_data$types == "tc")
idx_non_tc_obs = which(ntr_event_data$types == "non_tc")

for(i in 1:100){ 
  #Sample from multinomial
  type = rmultinom(1, size=44*5, prob = probs)
  
  #Number of each type
  n_tc = type[1]
  n_non_tc = type[2]
  
  #Pull out indicies of sampled rows
  rows = c(sample(idx_tc_obs,n_tc,replace=T),
           sample(idx_non_tc_obs,n_non_tc,replace=T))
  
  #Fit GPD
  gpd <- evm(obs[rows], th=min(obs[rows])) 
  gpd$rate <- n/44
  
  #Estimate retunr levels
  for(j in 1:100){
    rl_out = rl(gpd, M = gpd$rate*rp[j], alpha = 0.05, ci.fit = TRUE)[[1]]
    obs_rl_lower[i,j] = rl_out[2]
    obs_rl_est[i,j] = rl_out[1]
    obs_rl_upper[i,j] = rl_out[3]
  }
} 

write.csv(obs_rl_lower, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_lower.csv',sep=""))
write.csv(obs_rl_est, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_est.csv',sep=""))
write.csv(obs_rl_upper, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_upper.csv',sep=""))



  
  # ── wl marginal bootstrap ────────────────────────────────────────────────────

#Declustered water levels

#Read in CORA data
cora.file = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")
cora_wl_df = read.csv(cora.file)[,c(1,2)]
colnames(cora_wl_df) <- c("date", "wl")

cora_wl_df$date <- as.POSIXct(round((cora_wl_df$date - 719529) * 86400 / 3600) * 3600, 
                              origin = "1970-01-01", 
                              tz = "UTC")

#Add storm type
cora_wl_df <- HURDAT(cora_wl_df,
                     lat.loc = ba_lat_lon$lat,
                     lon.loc = -ba_lat_lon$lon,
                     rad = 350)

#Decluster ts
cora_wl_decl <- Decluster_SW(Data = cora_wl_df[,c(1,2)], Window_Width = 72)
cora_wl_decl_df <- data.frame(cora_wl_df$date,cora_wl_decl$Declustered)
colnames(cora_wl_decl_df) <- c("date", "wl")
cora_wl_decl_df$date <- as.Date(cora_wl_df$date)





# Event identification
ids <- order(cora_wl_decl_df$wl, decreasing = TRUE)[1:(5*44)]
types <- classify_event_types(ids, cora_wl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_wl_decl_df), ids + 1*24))
wl_event_data = list(ids = ids, types = types)

write.csv(wl_event_data,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_event_data_df.csv',sep=""), row.names = FALSE )

plot(cora_wl_df$date,cora_wl_df$wl)
points(cora_wl_df$date[wl_event_data$ids],cora_wl_df$wl[wl_event_data$ids], col=2)

#Write 3 WL components to an excel file
obs_wl_max_df = data.frame(cora_df[wl_event_data$ids,],"type" = wl_event_data$types)
write.csv(obs_wl_max_df,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_max_df.csv',sep="") )

events    <- df_sim$sample %in% c("ntr_tc", "ntr_non_tc")
event_idx <- which(events)

max_wl_id <- apply(wl$sim_wl_trim[event_idx, ], 1,
                   function(x) which.max(replace(x, is.na(x), -Inf)))

rc_idx <- cbind(event_idx, max_wl_id)

sim_wl_max_df <- data.frame(
  wl_max_wl   = wl$sim_wl_trim[rc_idx],
  wl_max_ntr  = wl$sim_ntr_trim[rc_idx],
  wl_max_pred = wl$sim_pred_trim[rc_idx],
  type        = df_sim$sample[event_idx]
)

write.csv(sim_wl_max_df, paste0('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/', name, '_sim_wl_max_df.csv'))


# Pre-compute quantities used repeatedly
sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
is_tc          <- wl_event_data$types == "tc"
rows_tc        <- which(df_sim$sample == "ntr_tc")
rows_non_tc    <- which(df_sim$sample == "ntr_non_tc")

wl_res <- list()
for (type in c("non_tc", "tc")) {
  
  obs      <- cora_wl_df$wl[wl_event_data$ids][if (type == "tc") is_tc else !is_tc]
  sim_rows <- if (type == "tc") rows_tc else rows_non_tc
  n        <- length(obs)
  rp       <- rp_axis(n)
  
  s               <- sort(obs)
  obs_rp10        <- approx(rp, s, xout = 10)$y
  obs_rp40         <- approx(rp, s, xout = 40)$y
  print("Here")
  gpd <- GPD_Fit(Data=s, Data_Full=cora_wl_df$wl, Thres=min(s), Method="Solari")
  gpd$Rate <- n/44
  
  #Extract ten year event
  obs_gpd_rp10  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((10*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp40  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((40*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp50  <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((50*gpd$Rate)^gpd$xi - 1)
  obs_gpd_rp100 <- gpd$Threshold + (gpd$sigma / gpd$xi) * ((100*gpd$Rate)^gpd$xi - 1)  
  
  obs_gpd_threshold  <- gpd$Threshold 
  obs_gpd_xi  <- gpd$xi
  obs_gpd_sigma <- gpd$sigma
  obs_ar2 <- AR2(c(obs_gpd_xi,obs_gpd_sigma,obs_gpd_threshold),obs)
  
  obs_perc = ((1:n)-0.5)/n
  
  bt <- bootstrap_rain(obs, sim_rows, sim_rain_sums=sim_wl_max, data=cora_wl_df$wl, perc=obs_perc)
  
  wl_res[[type]] <- list(
    obs = obs,
    obs_rp10           = obs_rp10,
    obs_rp40           = obs_rp40,
    obs_gpd_rp10       = obs_gpd_rp10, 
    obs_gpd_rp40       = obs_gpd_rp40, 
    obs_gpd_rp50       = obs_gpd_rp50,  
    obs_gpd_rp100      = obs_gpd_rp100,
    obs_gpd_threshold  = obs_gpd_threshold, 
    obs_gpd_xi         = obs_gpd_xi,
    obs_gpd_sigma      = obs_gpd_sigma, 
    obs_ar2            = obs_ar2, 
    rp10_025      = quantile(bt$rp10,0.025),
    rp10_05      = quantile(bt$rp10,0.05),
    rp10_25      = quantile(bt$rp10,0.25),
    rp10_50      = quantile(bt$rp10,0.5),
    rp10_75      = quantile(bt$rp10,0.75),
    rp10_95      = quantile(bt$rp10,0.95),
    rp10_975      = quantile(bt$rp10,0.975),
    rp40_025      = quantile(bt$rp40,0.025),
    rp40_05      = quantile(bt$rp40,0.05),
    rp40_25      = quantile(bt$rp40,0.25),
    rp40_50      = quantile(bt$rp40,0.5),
    rp40_75      = quantile(bt$rp40,0.75),
    rp40_95      = quantile(bt$rp40,0.95),
    rp40_975      = quantile(bt$rp40,0.975),
    gpd_rp10_025       = quantile(bt$gpd_rp10,0.025), 
    gpd_rp10_05       = quantile(bt$gpd_rp10,0.05), 
    gpd_rp10_25       = quantile(bt$gpd_rp10,0.25), 
    gpd_rp10_50       = quantile(bt$gpd_rp10,0.5),
    gpd_rp10_75       = quantile(bt$gpd_rp10,0.75),
    gpd_rp10_95       = quantile(bt$gpd_rp10,0.95), 
    gpd_rp10_975       = quantile(bt$gpd_rp10,0.975),
    gpd_rp40_025       = quantile(bt$gpd_rp40,0.025),
    gpd_rp40_05       = quantile(bt$gpd_rp40,0.05),
    gpd_rp40_25       = quantile(bt$gpd_rp40,0.25),
    gpd_rp40_50       = quantile(bt$gpd_rp40,0.50),
    gpd_rp40_75       = quantile(bt$gpd_rp40,0.75),
    gpd_rp40_95       = quantile(bt$gpd_rp40,0.95),
    gpd_rp40_975       = quantile(bt$gpd_rp40,0.975),
    gpd_rp50_025       = quantile(bt$gpd_rp50,0.025), 
    gpd_rp50_05       = quantile(bt$gpd_rp50,0.05), 
    gpd_rp50_25       = quantile(bt$gpd_rp50,0.25), 
    gpd_rp50_50       = quantile(bt$gpd_rp50,0.5),
    gpd_rp50_75       = quantile(bt$gpd_rp50,0.75),
    gpd_rp50_95       = quantile(bt$gpd_rp50,0.95),
    gpd_rp50_975       = quantile(bt$gpd_rp50,0.975), 
    gpd_rp100_025      = quantile(bt$gpd_rp100,0.025),
    gpd_rp100_05      = quantile(bt$gpd_rp100,0.05),
    gpd_rp100_25      = quantile(bt$gpd_rp100,0.25),
    gpd_rp100_50      = quantile(bt$gpd_rp100,0.50),
    gpd_rp100_75      = quantile(bt$gpd_rp100,0.75),
    gpd_rp100_95      = quantile(bt$gpd_rp100,0.95),
    gpd_rp100_975      = quantile(bt$gpd_rp100,0.975),
    gpd_threshold_025  = quantile(bt$gpd_threshold, 0.025),
    gpd_threshold_05  = quantile(bt$gpd_threshold, 0.05),
    gpd_threshold_25  = quantile(bt$gpd_threshold, 0.25),
    gpd_threshold_50  = quantile(bt$gpd_threshold, 0.50),
    gpd_threshold_75  = quantile(bt$gpd_threshold, 0.75),
    gpd_threshold_95  = quantile(bt$gpd_threshold, 0.95),
    gpd_threshold_975  = quantile(bt$gpd_threshold, 0.975),
    gpd_xi_025         = quantile(bt$gpd_xi,0.025),
    gpd_xi_05         = quantile(bt$gpd_xi,0.05),
    gpd_xi_25         = quantile(bt$gpd_xi,0.25),
    gpd_xi_50         = quantile(bt$gpd_xi,0.5),
    gpd_xi_75         = quantile(bt$gpd_xi,0.75),
    gpd_xi_95         = quantile(bt$gpd_xi,0.95),
    gpd_xi_975         = quantile(bt$gpd_xi,0.975),
    gpd_sigma_025      = quantile(bt$gpd_sigma, 0.025),
    gpd_sigma_05      = quantile(bt$gpd_sigma, 0.05),
    gpd_sigma_25      = quantile(bt$gpd_sigma, 0.25),
    gpd_sigma_50      = quantile(bt$gpd_sigma, 0.5),
    gpd_sigma_75      = quantile(bt$gpd_sigma, 0.75),
    gpd_sigma_95      = quantile(bt$gpd_sigma, 0.95),
    gpd_sigma_975      = quantile(bt$gpd_sigma, 0.975),
    ar2_025            = quantile(bt$ar2, 0.025),
    ar2_05            = quantile(bt$ar2, 0.05),
    ar2_25            = quantile(bt$ar2, 0.25),
    ar2_50            = quantile(bt$ar2, 0.5),
    ar2_75            = quantile(bt$ar2, 0.75),
    ar2_95            = quantile(bt$ar2, 0.95),
    ar2_975            = quantile(bt$ar2, 0.975),
    p_025            =  bt$p_025,
    p_05            =  bt$p_05,
    p_50            =  bt$p_50,
    p_95            =  bt$p_95,
    p_975            =  bt$p_975,
    p_val_reject = bt$p_val_reject
  )
}

res_wl<- data.frame(
  
  obs_rp10_tc           =  wl_res$tc$obs_rp10,
  obs_rp40_tc           =  wl_res$tc$obs_rp40,
  obs_gpd_rp10_tc       =  wl_res$tc$obs_gpd_rp10, 
  obs_gpd_rp40_tc       =  wl_res$tc$obs_gpd_rp40, 
  obs_gpd_rp50_tc       =  wl_res$tc$obs_gpd_rp50,  
  obs_gpd_rp100_tc      =  wl_res$tc$obs_gpd_rp100,
  obs_gpd_threshold_tc  =  wl_res$tc$obs_gpd_threshold, 
  obs_gpd_xi_tc         =  wl_res$tc$obs_gpd_xi,
  obs_gpd_sigma_tc      =  wl_res$tc$obs_gpd_sigma, 
  obs_ar2_tc            =  wl_res$tc$obs_ar2, 
  
  rp10_025_tc  = wl_res$tc$rp10_025,  
  rp10_05_tc  = wl_res$tc$rp10_05,          
  rp10_25_tc  =  wl_res$tc$rp10_25,           
  rp10_50_tc =  wl_res$tc$rp10_50,          
  rp10_75_tc  = wl_res$tc$rp10_75,      
  rp10_95_tc  =  wl_res$tc$rp10_95, 
  rp10_975_tc  = wl_res$tc$rp10_975,  
  
  rp40_025_tc  =  wl_res$tc$rp40_025, 
  rp40_05_tc  =   wl_res$tc$rp40_05, 
  rp40_25_tc  =  wl_res$tc$rp40_25,         
  rp40_50_tc  =  wl_res$tc$rp40_50,      
  rp40_75_tc  =   wl_res$tc$rp40_75,      
  rp40_95_tc  =  wl_res$tc$rp40_95,      
  rp40_975_tc  =  wl_res$tc$rp40_975, 
  
  gpd_rp10_025_tc  =  wl_res$tc$gpd_rp10_025,
  gpd_rp10_05_tc  =   wl_res$tc$gpd_rp10_05,        
  gpd_rp10_25_tc  =   wl_res$tc$gpd_rp10_25,        
  gpd_rp10_50_tc  =  wl_res$tc$gpd_rp10_50,      
  gpd_rp10_75_tc  =  wl_res$tc$gpd_rp10_75,       
  gpd_rp10_95_tc  =  wl_res$tc$gpd_rp10_95,
  gpd_rp10_975_tc  =  wl_res$tc$gpd_rp10_975,
  
  gpd_rp40_025_tc   =   wl_res$tc$gpd_rp40_025, 
  gpd_rp40_05_tc  =   wl_res$tc$gpd_rp40_05,       
  gpd_rp40_25_tc =   wl_res$tc$gpd_rp40_25,      
  gpd_rp40_50_tc   =   wl_res$tc$gpd_rp40_50,       
  gpd_rp40_75_tc   =   wl_res$tc$gpd_rp40_75,       
  gpd_rp40_95_tc   =   wl_res$tc$gpd_rp40_95, 
  gpd_rp40_975_tc   =   wl_res$tc$gpd_rp40_975, 
  
  gpd_rp50_025_tc  =   wl_res$tc$gpd_rp50_025, 
  gpd_rp50_05_tc  =   wl_res$tc$gpd_rp50_05,      
  gpd_rp50_25_tc  =   wl_res$tc$gpd_rp50_25,       
  gpd_rp50_50_tc  =   wl_res$tc$gpd_rp50_50,       
  gpd_rp50_75_tc  =   wl_res$tc$gpd_rp50_75,       
  gpd_rp50_95_tc  =   wl_res$tc$gpd_rp50_95, 
  gpd_rp50_975_tc  =   wl_res$tc$gpd_rp50_975, 
  
  gpd_rp100_025_tc  =   wl_res$tc$gpd_rp100_025,   
  gpd_rp100_05_tc  =   wl_res$tc$gpd_rp100_05,      
  gpd_rp100_25_tc  =   wl_res$tc$gpd_rp100_25,      
  gpd_rp100_50_tc  =   wl_res$tc$gpd_rp100_50,      
  gpd_rp100_75_tc  =   wl_res$tc$gpd_rp100_75,      
  gpd_rp100_95_tc  =   wl_res$tc$gpd_rp100_95,   
  gpd_rp100_975_tc  =   wl_res$tc$gpd_rp100_975,   
  
  gpd_threshold_025_tc  =   wl_res$tc$gpd_threshold_025,  
  gpd_threshold_05_tc  =   wl_res$tc$gpd_threshold_05,  
  gpd_threshold_25_tc  =   wl_res$tc$gpd_threshold_25, 
  gpd_threshold_50_tc  =   wl_res$tc$gpd_threshold_50,  
  gpd_threshold_75_tc  =   wl_res$tc$gpd_threshold_75,  
  gpd_threshold_95_tc  =   wl_res$tc$gpd_threshold_95,  
  gpd_threshold_975_tc  =   wl_res$tc$gpd_threshold_975,  
  
  gpd_xi_025_tc  =  wl_res$tc$gpd_xi_025, 
  gpd_xi_05_tc  =  wl_res$tc$gpd_xi_05,         
  gpd_xi_25_tc  =  wl_res$tc$gpd_xi_25,         
  gpd_xi_50_tc  =  wl_res$tc$gpd_xi_50,         
  gpd_xi_75_tc  =  wl_res$tc$gpd_xi_75,        
  gpd_xi_95_tc  =  wl_res$tc$gpd_xi_95, 
  gpd_xi_975_tc  =  wl_res$tc$gpd_xi_975, 
  
  gpd_sigma_025_tc  =   wl_res$tc$gpd_sigma_025,   
  gpd_sigma_05_tc  =   wl_res$tc$gpd_sigma_05,     
  gpd_sigma_25_tc  =   wl_res$tc$gpd_sigma_25,     
  gpd_sigma_50_tc  =   wl_res$tc$gpd_sigma_50,     
  gpd_sigma_75_tc  =   wl_res$tc$gpd_sigma_75,     
  gpd_sigma_95_tc  =   wl_res$tc$gpd_sigma_95,   
  gpd_sigma_975_tc  =   wl_res$tc$gpd_sigma_975,   
  
  ar2_025_tc    =  wl_res$tc$ar2_025,
  ar2_05_tc    =  wl_res$tc$ar2_05,            
  ar2_25_tc    =  wl_res$tc$ar2_25,            
  ar2_50_tc    =  wl_res$tc$ar2_50,            
  ar2_75_tc    =  wl_res$tc$ar2_75,            
  ar2_95_tc    =  wl_res$tc$ar2_95, 
  ar2_975_tc    =  wl_res$tc$ar2_975, 
  p_val_reject_tc  =   wl_res$tc$p_val_reject,
  
  obs_rp10_non_tc          =  wl_res$non_tc$obs_rp10,
  obs_rp40_non_tc          =  wl_res$non_tc$obs_rp40,
  obs_gpd_rp10_non_tc      =  wl_res$non_tc$obs_gpd_rp10, 
  obs_gpd_rp40_non_tc      =  wl_res$non_tc$obs_gpd_rp40, 
  obs_gpd_rp50_non_tc      =  wl_res$non_tc$obs_gpd_rp50,  
  obs_gpd_rp100_non_tc     =  wl_res$non_tc$obs_gpd_rp100,
  obs_gpd_threshold_non_tc  =  wl_res$non_tc$obs_gpd_threshold, 
  obs_gpd_xi_non_tc         =  wl_res$non_tc$obs_gpd_xi,
  obs_gpd_sigma_non_tc      =  wl_res$non_tc$obs_gpd_sigma, 
  obs_ar2_non_tc            =  wl_res$non_tc$obs_ar2, 
  
  rp10_025_non_tc  =  wl_res$non_tc$rp10_025,   
  rp10_05_non_tc  =  wl_res$non_tc$rp10_05,          
  rp10_25_non_tc  =  wl_res$non_tc$rp10_25,           
  rp10_50_non_tc  =  wl_res$non_tc$rp10_50,          
  rp10_75_non_tc  =  wl_res$non_tc$rp10_75,      
  rp10_95_non_tc  =  wl_res$non_tc$rp10_95,   
  rp10_975_non_tc  =  wl_res$non_tc$rp10_975,   
  
  rp40_025_non_tc  =   wl_res$non_tc$rp40_025,
  rp40_05_non_tc  =   wl_res$non_tc$rp40_05, 
  rp40_25_non_tc  =   wl_res$non_tc$rp40_25,         
  rp40_50_non_tc  =   wl_res$non_tc$rp40_50,      
  rp40_75_non_tc  =   wl_res$non_tc$rp40_75,      
  rp40_95_non_tc  =   wl_res$non_tc$rp40_95,
  rp40_975_non_tc  =   wl_res$non_tc$rp40_975,
  
  gpd_rp10_025_non_tc  =  wl_res$non_tc$gpd_rp10_025, 
  gpd_rp10_05_non_tc  =  wl_res$non_tc$gpd_rp10_05,        
  gpd_rp10_25_non_tc  =  wl_res$non_tc$gpd_rp10_25,        
  gpd_rp10_50_non_tc  =  wl_res$non_tc$gpd_rp10_50,      
  gpd_rp10_75_non_tc  =  wl_res$non_tc$gpd_rp10_75,       
  gpd_rp10_95_non_tc  =  wl_res$non_tc$gpd_rp10_95,        
  gpd_rp10_975_non_tc  =  wl_res$non_tc$gpd_rp10_975, 
  
  gpd_rp40_025_non_tc  =   wl_res$non_tc$gpd_rp40_025, 
  gpd_rp40_05_non_tc  =   wl_res$non_tc$gpd_rp40_05,       
  gpd_rp40_25_non_tc  =   wl_res$non_tc$gpd_rp40_25,      
  gpd_rp40_50_non_tc  =   wl_res$non_tc$gpd_rp40_50,       
  gpd_rp40_75_non_tc  =   wl_res$non_tc$gpd_rp40_75,       
  gpd_rp40_95_non_tc  =   wl_res$non_tc$gpd_rp40_95, 
  gpd_rp40_975_non_tc  =   wl_res$non_tc$gpd_rp40_975, 
  
  gpd_rp50_025_non_tc  =   wl_res$non_tc$gpd_rp50_025, 
  gpd_rp50_05_non_tc  =   wl_res$non_tc$gpd_rp50_05,      
  gpd_rp50_25_non_tc  =   wl_res$non_tc$gpd_rp50_25,       
  gpd_rp50_50_non_tc  =   wl_res$non_tc$gpd_rp50_50,       
  gpd_rp50_75_non_tc  =   wl_res$non_tc$gpd_rp50_75,       
  gpd_rp50_95_non_tc  =   wl_res$non_tc$gpd_rp50_95, 
  gpd_rp50_975_non_tc  =   wl_res$non_tc$gpd_rp50_975, 
  
  gpd_rp100_025_non_tc  =   wl_res$non_tc$gpd_rp100_025,  
  gpd_rp100_05_non_tc  =   wl_res$non_tc$gpd_rp100_05,      
  gpd_rp100_25_non_tc  =   wl_res$non_tc$gpd_rp100_25,      
  gpd_rp100_50_non_tc  =   wl_res$non_tc$gpd_rp100_50,      
  gpd_rp100_75_non_tc  =   wl_res$non_tc$gpd_rp100_75,      
  gpd_rp100_95_non_tc  =   wl_res$non_tc$gpd_rp100_95,  
  gpd_rp100_975_non_tc  =   wl_res$non_tc$gpd_rp100_975,  
  
  gpd_threshold_025_non_tc  =   wl_res$non_tc$gpd_threshold_025,  
  gpd_threshold_05_non_tc  =   wl_res$non_tc$gpd_threshold_05,  
  gpd_threshold_25_non_tc  =   wl_res$non_tc$gpd_threshold_25, 
  gpd_threshold_50_non_tc  =   wl_res$non_tc$gpd_threshold_50,  
  gpd_threshold_75_non_tc  =   wl_res$non_tc$gpd_threshold_75,  
  gpd_threshold_95_non_tc  =   wl_res$non_tc$gpd_threshold_95,  
  gpd_threshold_975_non_tc  =   wl_res$non_tc$gpd_threshold_975,  
  
  gpd_xi_025_non_tc  =  wl_res$non_tc$gpd_xi_025,         
  gpd_xi_05_non_tc  =  wl_res$non_tc$gpd_xi_05,         
  gpd_xi_25_non_tc  =  wl_res$non_tc$gpd_xi_25,         
  gpd_xi_50_non_tc  =  wl_res$non_tc$gpd_xi_50,         
  gpd_xi_75_non_tc  =  wl_res$non_tc$gpd_xi_75,        
  gpd_xi_95_non_tc  =  wl_res$non_tc$gpd_xi_95,         
  gpd_xi_975_non_tc  =  wl_res$non_tc$gpd_xi_975,         
  
  gpd_sigma_025_non_tc  =  wl_res$non_tc$gpd_sigma_025,
  gpd_sigma_05_non_tc  =   wl_res$non_tc$gpd_sigma_05,     
  gpd_sigma_25_non_tc  =   wl_res$non_tc$gpd_sigma_25,     
  gpd_sigma_50_non_tc  =   wl_res$non_tc$gpd_sigma_50,     
  gpd_sigma_75_non_tc  =   wl_res$non_tc$gpd_sigma_75,     
  gpd_sigma_95_non_tc  =  wl_res$non_tc$gpd_sigma_95,
  gpd_sigma_975_non_tc  =  wl_res$non_tc$gpd_sigma_975,
  
  ar2_025_non_tc   =  wl_res$non_tc$ar2_025,  
  ar2_05_non_tc   =  wl_res$non_tc$ar2_05,            
  ar2_25_non_tc   =  wl_res$non_tc$ar2_25,            
  ar2_50_non_tc   =  wl_res$non_tc$ar2_50,            
  ar2_75_non_tc   =  wl_res$non_tc$ar2_75,            
  ar2_95_non_tc   =  wl_res$non_tc$ar2_95,
  ar2_975_non_tc   =  wl_res$non_tc$ar2_975,  
  
  p_val_reject_non_tc = wl_res$non_tc$p_val_reject)

res_wl_df_tc= data.frame(obs_tc             =  wl_res$tc$obs,
                         p025_tc            =  wl_res$tc$p_025,
                         p05_tc             =  wl_res$tc$p_05,
                         p50_tc             =  wl_res$tc$p_50,
                         p95_tc             =  wl_res$tc$p_95,
                         p975_tc            =  wl_res$tc$p_975)

res_wl_df_non_tc  = data.frame(obs_non_tc            =  wl_res$non_tc$obs,
                               p025_non_tc            =  wl_res$non_tc$p_025,
                               p05_non_tc            =  wl_res$non_tc$p_05,
                               p50_non_tc            =  wl_res$non_tc$p_50,
                               p95_non_tc            =  wl_res$non_tc$p_95,
                               p975_non_tc           =  wl_res$non_tc$p_975)


write.csv(res_wl, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_res_aorc_subset.csv',sep=""))
write.csv(res_wl_df_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_res_df_tc_aorc_subset.csv',sep=""))
write.csv(res_wl_df_non_tc, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_res_df_non_tc_aorc_subset.csv',sep=""))


obs = cora_wl_df$wl[wl_event_data$ids]
n = length(obs)

#Probabilities for multinomial distribution
total_events = 5*44
p_ntr_tc = length(which(wl_event_data$types == "tc"))/total_events
p_ntr_non_tc = length(which(wl_event_data$types == "non_tc"))/total_events
probs = c(p_ntr_tc, p_ntr_non_tc)

obs_rl_lower = array(NA, dim=c(100,100))
obs_rl_est = array(NA, dim=c(100,100))
obs_rl_upper = array(NA, dim=c(100,100))

rp= 1:100

#Pull out indicies of TC and non-tc conditioned 
idx_tc_obs  = which(wl_event_data$types == "tc")
idx_non_tc_obs = which(wl_event_data$types == "non_tc")

for(i in 1:100){ 
  #Sample from multinomial
  type = rmultinom(1, size=44*5, prob = probs)
  
  #Number of each type
  n_tc = type[1]
  n_non_tc = type[2]
  
  #Pull out indicies of sampled rows
  rows = c(sample(idx_tc_obs,n_tc, replace=T),
           sample(idx_non_tc_obs,n_non_tc, replace=T))
  
  #Fit GPD
  gpd <- evm(obs[rows], th=min(obs[rows])) 
  gpd$rate <- n/44
  
  #Estimate retunr levels
  for(j in 1:100){
    rl_out = rl(gpd, M = gpd$rate*rp[j], alpha = 0.05, ci.fit = TRUE)[[1]]
    obs_rl_lower[i,j] = rl_out[2]
    obs_rl_est[i,j] = rl_out[1]
    obs_rl_upper[i,j] = rl_out[3]
  }
} 

write.csv(obs_rl_lower, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_lower.csv',sep=""))
write.csv(obs_rl_est, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_est.csv',sep=""))
write.csv(obs_rl_upper, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_upper.csv',sep=""))

#Simulations
sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)

#Observations
obs      <- cora_wl_df$wl[wl_event_data$ids]
n        <- length(obs)

#Probabilities for multinomial distribution
total_events = 5*44
p_ntr_tc = length(which(wl_event_data$types == "tc"))/total_events
p_ntr_non_tc = length(which(wl_event_data$types == "non_tc"))/total_events
probs = c(p_ntr_tc, p_ntr_non_tc)

sim_rl_lower = array(NA, dim=c(100,100))
sim_rl_est = array(NA, dim=c(100,100))
sim_rl_upper = array(NA, dim=c(100,100))

rp= 1:100

#Pull out indicies of Tc and non-TC events
idx_tc_sim = which(df_sim$sample == "ntr_tc")
idx_non_tc_sim = which(df_sim$sample == "ntr_non_tc")     

for(i in 1:100){ 
  #Sample from multinomial
  type = rmultinom(1, size=44*5, prob = probs)
  
  #Pull out number of sampled Tc and non-TC events
  n_tc = type[1]
  n_non_tc = type[2]       
  
  #Obtain the sample
  rows = c(sample(idx_tc_sim,n_tc,replace=T),
           sample(idx_non_tc_sim,n_non_tc,replace=T))   
  
  #Fit GPD
  gpd <- evm(sim_wl_max[rows], th=min(sim_wl_max[rows])) 
  gpd$rate <- n/44
  
  #Estimate return levels
  for(j in 1:100){
    rl_out = rl(gpd, M = gpd$rate*rp[j], alpha = 0.05, ci.fit = TRUE)[[1]]
    sim_rl_lower[i,j] = rl_out[2]
    sim_rl_est[i,j] = rl_out[1]
    sim_rl_upper[i,j] = rl_out[3]
  }
} 

write.csv(sim_rl_lower, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_lower.csv',sep=""))
write.csv(sim_rl_est, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_est.csv',sep=""))
write.csv(sim_rl_upper, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_upper.csv',sep=""))


elapsed <- proc.time() - start
print(elapsed)

write.csv(wl$sim_wl_trim, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_wl_trim.csv',sep=""))
write.csv(wl$sim_ntr_trim, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_ntr_trim.csv',sep=""))
write.csv(wl$sim_pred_trim, paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_pred_trim.csv',sep=""))



vol = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_aorc_subset.csv',sep=""))[,-1]
vol_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_df_tc_aorc_subset.csv',sep=""))[,-1]
vol_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_vol_df_non_tc_aorc_subset.csv',sep=""))[,-1]

plot(vol$obs_gpd_sigma_tc,vol$obs_gpd_xi_tc,xlim=c(0,100),ylim=c(-1,1))
points(vol$gpd_sigma_50_tc,vol$gpd_xi_50_tc,col=2)

diff_10  = vol$obs_rp10_tc  - vol$rp10_50_tc
diff_40  = vol$obs_rp40_tc  - vol$rp40_50_tc
diff_gpd_10  = vol$obs_gpd_rp10_tc  - vol$gpd_rp10_50_tc
diff_gpd_50  = vol$obs_gpd_rp50_tc  - vol$gpd_rp50_50_tc
diff_gpd_100 = vol$obs_gpd_rp100_tc - vol$gpd_rp100_50_tc

diff_10_non_tc   = vol$obs_rp10_non_tc  - vol$rp10_50_non_tc 
diff_40_non_tc   = vol$obs_rp40_non_tc   - vol$rp40_50_non_tc 
diff_gpd_10_non_tc   = vol$obs_gpd_rp10_non_tc   - vol$gpd_rp10_50_non_tc 
diff_gpd_50_non_tc   = vol$obs_gpd_rp50_non_tc   - vol$gpd_rp50_50_non_tc 
diff_gpd_100_non_tc  = vol$obs_gpd_rp100_non_tc  - vol$gpd_rp100_50_non_tc 

data_list <- list(diff_10, diff_10_non_tc, diff_40, diff_40_non_tc, diff_gpd_10,diff_gpd_10_non_tc, diff_gpd_50,diff_gpd_50_non_tc, diff_gpd_100, diff_gpd_100_non_tc)
labels    <- c("10-year (empirical)","40-year (empirical)","10-year GPD","50-year GPD","100-year GPD")
box_cols  <- c("red","blue","red","blue","red","blue","red","blue","red","blue")


png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Obs_minus_sim_rain_vol_boxplots.png",width=10,height=7.5,res=300,units='in')
par(mar=c(4.2,4.5,0.5,0.5))
# ── Function: compute box stats for one vector ────────────────
box_stats <- function(x) {
  q1     <- quantile(x, 0.25)
  med    <- median(x)
  q3     <- quantile(x, 0.75)
  iqr    <- q3 - q1
  w_low  <- min(x[x >= q1 - 1.5 * iqr])   # lower whisker end
  w_high <- max(x[x <= q3 + 1.5 * iqr])   # upper whisker end
  out    <- x[x < w_low | x > w_high]      # outliers
  list(q1=q1, med=med, q3=q3, w_low=w_low, w_high=w_high, out=out)
}

stats <- lapply(data_list, box_stats)

# ── Set up blank plot ─────────────────────────────────────────
all_vals <- unlist(data_list)
plot(0, type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(min(all_vals) - 5, max(all_vals) + 5),
     xlab="",ylab="",
     xaxt = "n")                            # suppress default x-axis
abline(h=0)
mtext("Event", side=1,line=2.5)
mtext("Observed - Simulated (mm)", side = 2, line = 2.5)

axis(1, at = seq(1,10,2), labels = labels)          # custom x-axis labels
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray80")

# ── Draw each box-and-whisker manually ───────────────────────
box_width <- 0.35/2                           # half-width of box
cap_width <- 0.15/2                           # half-width of whisker caps

for (i in 1:10) {
  s  <- stats[[i]]
  cx <- c(0.75,1.25,2.75,3.25,4.75,5.25,6.75,7.25,8.75,9.25)[i]                                 # centre x position
  
  # 1. Box (Q1 to Q3)
  rect(xleft   = cx - box_width,
       ybottom = s$q1,
       xright  = cx + box_width,
       ytop    = s$q3,
       col     = box_cols[i],
       border  = "black",
       lwd     = 1.5)
  
  # 2. Median line
  segments(x0  = cx - box_width, y0 = s$med,
           x1  = cx + box_width, y1 = s$med,
           lwd = 2.5, col = "black")
  
  # 3. Whisker lines (centre of box to whisker ends)
  segments(x0 = cx, y0 = s$q1,    x1 = cx, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx, y0 = s$q3,    x1 = cx, y1 = s$w_high, lwd = 1.5)
  
  # 4. Whisker caps (horizontal bars at each end)
  segments(x0 = cx - cap_width, y0 = s$w_low,
           x1 = cx + cap_width, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx - cap_width, y0 = s$w_high,
           x1 = cx + cap_width, y1 = s$w_high, lwd = 1.5)
  
  # 5. Outlier points
  if (length(s$out) > 0) {
    points(rep(cx, length(s$out)), s$out,
           pch = 16, cex = 0.85, col = box_cols[i])
  }
}

dev.off()



png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Obs_minus_sim_rain_vol_boxplots.png",width=10,height=7.5,res=300,units='in')
par(mar=c(4.2,4.5,0.5,0.5))
# ── Function: compute box stats for one vector ────────────────
box_stats <- function(x) {
  q1     <- quantile(x, 0.25)
  med    <- median(x)
  q3     <- quantile(x, 0.75)
  iqr    <- q3 - q1
  w_low  <- min(x[x >= q1 - 1.5 * iqr])   # lower whisker end
  w_high <- max(x[x <= q3 + 1.5 * iqr])   # upper whisker end
  out    <- x[x < w_low | x > w_high]      # outliers
  list(q1=q1, med=med, q3=q3, w_low=w_low, w_high=w_high, out=out)
}

stats <- lapply(data_list, box_stats)

# ── Set up blank plot ─────────────────────────────────────────
all_vals <- unlist(data_list)
plot(0, type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(min(all_vals) - 5, max(all_vals) + 5),
     xlab="",ylab="",
     xaxt = "n")                            # suppress default x-axis
abline(h=0)
mtext("Event", side=1,line=2.5)
mtext("Observed - Simulated (mm)", side = 2, line = 2.5)

axis(1, at = seq(1,10,2), labels = labels)          # custom x-axis labels
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray80")

# ── Draw each box-and-whisker manually ───────────────────────
box_width <- 0.35/2                           # half-width of box
cap_width <- 0.15/2                           # half-width of whisker caps

for (i in 1:10) {
  s  <- stats[[i]]
  cx <- c(0.75,1.25,2.75,3.25,4.75,5.25,6.75,7.25,8.75,9.25)[i]                                 # centre x position
  
  # 1. Box (Q1 to Q3)
  rect(xleft   = cx - box_width,
       ybottom = s$q1,
       xright  = cx + box_width,
       ytop    = s$q3,
       col     = box_cols[i],
       border  = "black",
       lwd     = 1.5)
  
  # 2. Median line
  segments(x0  = cx - box_width, y0 = s$med,
           x1  = cx + box_width, y1 = s$med,
           lwd = 2.5, col = "black")
  
  # 3. Whisker lines (centre of box to whisker ends)
  segments(x0 = cx, y0 = s$q1,    x1 = cx, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx, y0 = s$q3,    x1 = cx, y1 = s$w_high, lwd = 1.5)
  
  # 4. Whisker caps (horizontal bars at each end)
  segments(x0 = cx - cap_width, y0 = s$w_low,
           x1 = cx + cap_width, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx - cap_width, y0 = s$w_high,
           x1 = cx + cap_width, y1 = s$w_high, lwd = 1.5)
  
  # 5. Outlier points
  if (length(s$out) > 0) {
    points(rep(cx, length(s$out)), s$out,
           pch = 16, cex = 0.85, col = box_cols[i])
  }
}

dev.off()



wl = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_res_aorc_subset.csv',sep=""))[,-1]

diff_10  = wl$obs_rp10_tc  - wl$rp10_50_tc
diff_40  = wl$obs_rp40_tc  - wl$rp40_50_tc
diff_gpd_10  = wl$obs_gpd_rp10_tc  - wl$gpd_rp10_50_tc
diff_gpd_50  = wl$obs_gpd_rp50_tc  - wl$gpd_rp50_50_tc
diff_gpd_100 = wl$obs_gpd_rp100_tc - wl$gpd_rp100_50_tc

diff_10_non_tc   = wl$obs_rp10_non_tc  -  wl$rp10_50_non_tc 
diff_40_non_tc   = wl$obs_rp40_non_tc   - wl$rp40_50_non_tc 
diff_gpd_10_non_tc   = wl$obs_gpd_rp10_non_tc   - wl$gpd_rp10_50_non_tc 
diff_gpd_50_non_tc   = wl$obs_gpd_rp50_non_tc   - wl$gpd_rp50_50_non_tc 
diff_gpd_100_non_tc  = wl$obs_gpd_rp100_non_tc  - wl$gpd_rp100_50_non_tc 

data_list <- list(diff_10, diff_10_non_tc, diff_40, diff_40_non_tc, diff_gpd_10,diff_gpd_10_non_tc, diff_gpd_50,diff_gpd_50_non_tc, diff_gpd_100, diff_gpd_100_non_tc)
labels    <- c("10-year (empirical)","40-year (empirical)","10-year GPD","50-year GPD","100-year GPD")
box_cols  <- c("red","blue","red","blue","red","blue","red","blue","red","blue")

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Obs_minus_sim_wl_boxplots.png",width=10,height=7.5,res=300,units='in')
par(mar=c(4.2,4.5,0.5,0.5))
# ── Function: compute box stats for one vector ────────────────
box_stats <- function(x) {
  q1     <- quantile(x, 0.25)
  med    <- median(x)
  q3     <- quantile(x, 0.75)
  iqr    <- q3 - q1
  w_low  <- min(x[x >= q1 - 1.5 * iqr])   # lower whisker end
  w_high <- max(x[x <= q3 + 1.5 * iqr])   # upper whisker end
  out    <- x[x < w_low | x > w_high]      # outliers
  list(q1=q1, med=med, q3=q3, w_low=w_low, w_high=w_high, out=out)
}

stats <- lapply(data_list, box_stats)

# ── Set up blank plot ─────────────────────────────────────────
all_vals <- unlist(data_list)
plot(0, type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(min(all_vals) - 0.1, max(all_vals) + 0.1),
     xlab="",ylab="",
     xaxt = "n")                            # suppress default x-axis
abline(h=0)
mtext("Event", side=1,line=2.5)
mtext("Observed - Simulated (m)", side = 2, line = 2.5)

axis(1, at = seq(1,10,2), labels = labels)          # custom x-axis labels
abline(h = seq(-0.6,0.6,0.2), lty = "dotted", col = "gray80")

# ── Draw each box-and-whisker manually ───────────────────────
box_width <- 0.35/2                           # half-width of box
cap_width <- 0.15/2                           # half-width of whisker caps

for (i in 1:10) {
  s  <- stats[[i]]
  cx <- c(0.75,1.25,2.75,3.25,4.75,5.25,6.75,7.25,8.75,9.25)[i]                                 # centre x position
  
  # 1. Box (Q1 to Q3)
  rect(xleft   = cx - box_width,
       ybottom = s$q1,
       xright  = cx + box_width,
       ytop    = s$q3,
       col     = box_cols[i],
       border  = "black",
       lwd     = 1.5)
  
  # 2. Median line
  segments(x0  = cx - box_width, y0 = s$med,
           x1  = cx + box_width, y1 = s$med,
           lwd = 2.5, col = "black")
  
  # 3. Whisker lines (centre of box to whisker ends)
  segments(x0 = cx, y0 = s$q1,    x1 = cx, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx, y0 = s$q3,    x1 = cx, y1 = s$w_high, lwd = 1.5)
  
  # 4. Whisker caps (horizontal bars at each end)
  segments(x0 = cx - cap_width, y0 = s$w_low,
           x1 = cx + cap_width, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx - cap_width, y0 = s$w_high,
           x1 = cx + cap_width, y1 = s$w_high, lwd = 1.5)
  
  # 5. Outlier points
  if (length(s$out) > 0) {
    points(rep(cx, length(s$out)), s$out,
           pch = 16, cex = 0.85, col = box_cols[i])
  }
}

dev.off()



plot(vol$obs_gpd_sigma_tc,vol$obs_gpd_xi_tc,xlim=c(0,100),ylim=c(-1,1))
points(vol$gpd_sigma_50_tc,vol$gpd_xi_50_tc,col=2)



ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_wl_res_aorc_subset.csv',sep=""))[,-1]



rain_vol = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_aorc_subset.csv',sep=""))[,-1]
rain_vol_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_tc_aorc_subset.csv',sep=""))[,-1]
rain_vol_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_non_tc_aorc_subset.csv',sep=""))[,-1]

diff_10  = rain_vol$obs_rp10_tc  - rain_vol$rp10_50_tc
diff_40  = rain_vol$obs_rp40_tc  - rain_vol$rp40_50_tc
diff_gpd_10  = rain_vol$obs_gpd_rp10_tc  - rain_vol$gpd_rp10_50_tc
diff_gpd_50  = rain_vol$obs_gpd_rp50_tc  - rain_vol$gpd_rp50_50_tc
diff_gpd_100 = rain_vol$obs_gpd_rp100_tc - rain_vol$gpd_rp100_50_tc

diff_10_non_tc   = rain_vol$obs_rp10_non_tc  - rain_vol$rp10_50_non_tc 
diff_40_non_tc   = rain_vol$obs_rp40_non_tc   - rain_vol$rp40_50_non_tc 
diff_gpd_10_non_tc   = rain_vol$obs_gpd_rp10_non_tc   - rain_vol$gpd_rp10_50_non_tc 
diff_gpd_50_non_tc   = rain_vol$obs_gpd_rp50_non_tc   - rain_vol$gpd_rp50_50_non_tc 
diff_gpd_100_non_tc  = rain_vol$obs_gpd_rp100_non_tc  - rain_vol$gpd_rp100_50_non_tc 

data_list <- list(diff_10, diff_10_non_tc, diff_40, diff_40_non_tc, diff_gpd_10,diff_gpd_10_non_tc, diff_gpd_50,diff_gpd_50_non_tc, diff_gpd_100, diff_gpd_100_non_tc)
labels    <- c("10-year (empirical)","40-year (empirical)","10-year GPD","50-year GPD","100-year GPD")
box_cols  <- c("red","blue","red","blue","red","blue","red","blue","red","blue")

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/rain_res_boxplots.png",width=10,height=7.5,res=300,units='in')

# ── Function: compute box stats for one vector ────────────────
box_stats <- function(x) {
  q1     <- quantile(x, 0.25)
  med    <- median(x)
  q3     <- quantile(x, 0.75)
  iqr    <- q3 - q1
  w_low  <- min(x[x >= q1 - 1.5 * iqr])   # lower whisker end
  w_high <- max(x[x <= q3 + 1.5 * iqr])   # upper whisker end
  out    <- x[x < w_low | x > w_high]      # outliers
  list(q1=q1, med=med, q3=q3, w_low=w_low, w_high=w_high, out=out)
}

stats <- lapply(data_list, box_stats)

# ── Set up blank plot ─────────────────────────────────────────
all_vals <- unlist(data_list)
plot(0, type = "n",
     xlim = c(0.5, 9.5),
     ylim = c(min(all_vals) - 5, max(all_vals) + 5),
     xlab="",ylab="",
     xaxt = "n")                            # suppress default x-axis
abline(h=0)
mtext("Event", side=1,line=2.5)
mtext("Observed - Simulated (mm)", side = 2, line = 2.5)

axis(1, at = seq(1,10,2), labels = labels)          # custom x-axis labels
grid(nx = NA, ny = NULL, lty = "dotted", col = "gray80")

# ── Draw each box-and-whisker manually ───────────────────────
box_width <- 0.35/2                           # half-width of box
cap_width <- 0.15/2                           # half-width of whisker caps

for (i in 1:10) {
  s  <- stats[[i]]
  cx <- c(0.75,1.25,2.75,3.25,4.75,5.25,6.75,7.25,8.75,9.25)[i]                                 # centre x position
  
  # 1. Box (Q1 to Q3)
  rect(xleft   = cx - box_width,
       ybottom = s$q1,
       xright  = cx + box_width,
       ytop    = s$q3,
       col     = box_cols[i],
       border  = "black",
       lwd     = 1.5)
  
  # 2. Median line
  segments(x0  = cx - box_width, y0 = s$med,
           x1  = cx + box_width, y1 = s$med,
           lwd = 2.5, col = "black")
  
  # 3. Whisker lines (centre of box to whisker ends)
  segments(x0 = cx, y0 = s$q1,    x1 = cx, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx, y0 = s$q3,    x1 = cx, y1 = s$w_high, lwd = 1.5)
  
  # 4. Whisker caps (horizontal bars at each end)
  segments(x0 = cx - cap_width, y0 = s$w_low,
           x1 = cx + cap_width, y1 = s$w_low,  lwd = 1.5)
  segments(x0 = cx - cap_width, y0 = s$w_high,
           x1 = cx + cap_width, y1 = s$w_high, lwd = 1.5)
  
  # 5. Outlier points
  if (length(s$out) > 0) {
    points(rep(cx, length(s$out)), s$out,
           pch = 1, cex = 1.2, col = box_cols[i])
  }
}

dev.off()


res_rain_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_tc_aorc_subset_1.csv',sep=""))
plot(0,type='n',xlim=c(10,300),ylim=c(10,300))
for(i in 1:10){
  points(sort(res_rain_tc$obs_tc[which(res_rain_tc$aorc==i)]),
         res_rain_tc$p50_tc[which(res_rain_tc$aorc==i)])
}
lines(seq(-10,350,0.1),seq(-10,350,0.1))



res_rain_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_tc_aorc_subset_1.csv',sep=""))
plot(0,type='n',xlim=c(10,300),ylim=c(10,300))
for(i in 1:10){
  points(sort(res_rain_tc$obs_tc[which(res_rain_tc$aorc==i)]),
         res_rain_tc$p50_tc[which(res_rain_tc$aorc==i)])
}
lines(seq(-10,350,0.1),seq(-10,350,0.1))

res_rain_non_tc = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_rain_res_df_non_tc_aorc_subset_1.csv',sep=""))
plot(0,type='n',xlim=c(10,300),ylim=c(10,300))
for(i in 1:10){
  points(sort(res_rain_non_tc$obs_non_tc[which(res_rain_non_tc$aorc==i)]),
         res_rain_non_tc$p50_non_tc[which(res_rain_non_tc$aorc==i)])
}
lines(seq(-10,350,0.1),seq(-10,350,0.1))


#Check ntr
sim_wl_max     <- apply(wl$sim_wl_trim, 1, max)
sim_ntr_max     <- apply(wl$sim_ntr_trim, 1, max)
df_sim$ntr


sum(ifelse(round(sim_ntr_max,3) == round(max(df_sim$ntr),3),1,0))


##GPD ntr obs vs sim

obs_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_lower.csv',sep=""))[,-1]
obs_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_est.csv',sep=""))[,-1]
obs_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_gpd_combined_upper.csv',sep=""))[,-1]

sim_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_lower.csv',sep=""))[,-1]
sim_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_est.csv',sep=""))[,-1]
sim_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_ntr_gpd_combined_upper.csv',sep=""))[,-1]

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/ntr_gpd_return_level.png",width=10,height=7.5,res=300,units='in')

plot(0, type='n', xlim=c(1,100), ylim=c(0.5,3),
     xaxt='n',yaxt='n',xlab="", ylab="")

mtext("Return Period (years)", side=1,line=2.75)
mtext("NTR (m)", side = 2, line = 2.75)

axis(1, at =  c(1,20,40,60,80,100), labels =   c(1,20,40,60,80,100))      
axis(1, at = seq(10,90,20), labels = FALSE, tck= -0.005)  

axis(2, at = seq(0.5,3,0.5), labels = seq(0.5,3,0.5))    
axis(2, at = seq(0.75,2.75,0.25), labels = FALSE, tck=-0.005)  

# ── Obs envelope ──────────────────────────────────────────────
obs_lower_env <- apply(obs_est, 2, quantile, 0.025)
obs_upper_env <- apply(obs_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(obs_upper_env, rev(obs_lower_env)),
        col = adjustcolor("black", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(obs_est, 2, mean), col = "black", lwd = 2)

# ── Sim envelope ──────────────────────────────────────────────
sim_lower_env <- apply(sim_est, 2, quantile, 0.025)
sim_upper_env <- apply(sim_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(sim_upper_env, rev(sim_lower_env)),
        col = adjustcolor("red", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(sim_est, 2, mean), col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed", "Simulated"),
       lty    = 1, lwd = 2, col = c("black","red"),
       bty    = "n")

dev.off()

#WL
obs_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_lower.csv',sep=""))[,-1]
obs_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_est.csv',sep=""))[,-1]
obs_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_gpd_combined_upper.csv',sep=""))[,-1]

sim_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_lower.csv',sep=""))[,-1]
sim_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_est.csv',sep=""))[,-1]
sim_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_upper.csv',sep=""))[,-1]

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/gpd_wl_return_level.png",width=10,height=7.5,res=300,units='in')

plot(0, type='n', xlim=c(1,100), ylim=c(1,3.5),
     xaxt='n',yaxt='n',xlab="", ylab="")

mtext("Return Period (years)", side=1,line=2.75)
mtext("WL (m)", side = 2, line = 2.75)

axis(1, at =  c(1,20,40,60,80,100), labels =   c(1,20,40,60,80,100))      
axis(1, at = seq(10,90,20), labels = FALSE, tck= -0.005)  

axis(2, at = seq(1,3.5,0.5), labels = seq(1,3.5,0.5))    
axis(2, at = seq(1.25,3.25,0.25), labels = FALSE, tck=-0.005)  

# ── Obs envelope ──────────────────────────────────────────────
obs_lower_env <- apply(obs_est, 2, quantile, 0.025)
obs_upper_env <- apply(obs_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(obs_upper_env, rev(obs_lower_env)),
        col = adjustcolor("black", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(obs_est, 2, mean), col = "black", lwd = 2)

# ── Sim envelope ──────────────────────────────────────────────
sim_lower_env <- apply(sim_est, 2, quantile, 0.025)
sim_upper_env <- apply(sim_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(sim_upper_env, rev(sim_lower_env)),
        col = adjustcolor("red", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(sim_est, 2, mean), col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed", "Simulated"),
       lty    = 1, lwd = 2, col = c("black","red"),
       bty    = "n")


dev.off()


##Water level analysis

#Find location of ntrmaximum
wl_tc_max_id = max.col(wl$sim_wl_trim[which(df_sim$sample=='ntr_tc'),]) 
ntr_tc_max_id = max.col(wl$sim_ntr_trim[which(df_sim$sample=='ntr_tc'),])
length(which(wl_tc_max_id==37))/length(which(df_sim$sample=='ntr_tc')) #17%

wl_non_tc_max_id = max.col(wl$sim_wl_trim[which(df_sim$sample=='ntr_non_tc'),]) 
ntr_non_tc_max_id = max.col(wl$sim_ntr_trim[which(df_sim$sample=='ntr_non_tc'),])
length(which(wl_non_tc_max_id==37))/length(which(df_sim$sample=='ntr_non_tc')) #11%

#See if it matches maximum
plot(df_sim$ntr_factor[which(df_sim$sample=='ntr_non_tc')],apply(wl$sim_wl_trim,1,max)[which(df_sim$sample=='ntr_non_tc')], col=ifelse(wl_non_tc_max_id==37,"Red","Black"))
points(df_sim$ntr_factor[which(df_sim$sample=='ntr_tc')],apply(wl$sim_wl_trim,1,max)[which(df_sim$sample=='ntr_tc')], col=ifelse(wl_tc_max_id==37,"Red","Black"))

#Compare correlations in simulated and observed samples rainfall - NTR
cor_r = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_con_rain_aorc_subset.csv',sep=""))
cor_ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_con_ntr_aorc_subset.csv',sep=""))

cor_r = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_24hr_rain_aorc_subset.csv',sep=""))
cor_ntr = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_ntr_con_ntr_aorc_subset.csv',sep=""))

cor_obs = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Mullica_cor_by_acc_time_',50,'.csv',sep=""))
n_files = seq(100,2300,50)[1:10]
for(i in 2:length(n_files)){
  print(i)
  new  = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Mullica_cor_by_acc_time_',n_files[i],'.csv',sep=""))
  cor_obs = rbind(cor_obs, new)
}


z = which(cor_obs$rain_dur_seq==24)

par(mfrow=c(2,2))
plot(cor_obs$cor_r_tc[z][1:500],cor_r$cor_r_tc[1:500],xlim=c(0,0.5),ylim=c(0,0.5))#,xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
lines(seq(0,1,0.1),seq(0,1,0.1))
plot(cor_obs$cor_r_non_tc[z][1:500],cor_r$cor_r_non_tc[1:500],xlim=c(-0.1,0.1),ylim=c(-0.1,0.1))
lines(seq(0,1,0.1),seq(0,1,0.1))

plot(cor_obs$cor_ntr_tc[z][1:500],cor_ntr$cor_ntr_tc[1:500],xlim=c(0.1,0.4),ylim=c(0.1,0.4))
lines(seq(0,1,0.1),seq(0,1,0.1))
plot(cor_obs$cor_ntr_non_tc[z][1:500],cor_ntr$cor_ntr_non_tc[1:500],xlim=c(0.1,0.4),ylim=c(0.1,0.4))
lines(seq(0,1,0.1),seq(0,1,0.1))

#Compare correlations in simulated and observed samples rainfall -WL
cor_r = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_con_rain_aorc_subset.csv',sep=""))
cor_wl = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_rain_wl_con_wl_aorc_subset.csv',sep=""))

cor_r = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_wl_con_24hr_rain_aorc_subset.csv',sep=""))
cor_wl = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_res_cor_24hr_rain_wl_con_wl_aorc_subset.csv',sep=""))

cor_obs = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/Mullica_anv_spatial_24_hr_acc_time_10_yr_r_wl_',100,'.csv',sep=""))
n_files = c(seq(200,2000,100),2300)
for(i in 1:length(n_files)){
  print(i)
  new  = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/con_wl/data/Mullica_anv_spatial_24_hr_acc_time_10_yr_r_wl_',n_files[i],'.csv',sep=""))
  cor_obs = rbind(cor_obs, new)
}


par(mfrow=c(2,2))
plot(cor_obs$cor_r_tc[1:500],cor_r$cor_r_tc[1:500],xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
plot(cor_obs$cor_r_non_tc[1:500],cor_r$cor_r_non_tc[1:500],xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))

plot(cor_obs$cor_wl_tc[1:500],cor_wl$cor_wl_tc[1:500],xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))
plot(cor_obs$cor_wl_non_tc[1:500],cor_wl$cor_wl_non_tc[1:500],xlim=c(-0.5,0.5),ylim=c(-0.5,0.5))





##DDF analysis QQ-plots
df_tc = read.csv('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/Mullica_ddf_tc_aorc_subset.csv')
df_non_tc = read.csv('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/Mullica_ddf_non_tc_aorc_subset.csv')

event = data.frame(
  dur = rep(c(12,24,36), 5),
  rp  = rep(c(5,10,25,50,100), each = 3)
)

# Layout:
# extra top row for duration labels
# extra left column for RP labels
lay <- matrix(c(
  1,  2,  3,  4,
  5,  6,  7,  8,
  9, 10, 11, 12,
  13, 14, 15, 16,
  17, 18, 19, 20,
  21, 22, 23, 24
), nrow = 6, byrow = TRUE)

layout(lay, widths = c(1,4,4,4), heights = c(1,4,4,4,4,4))

par(mar=c(0,0,0,0))

# Empty top-left corner
plot.new()

# Top duration labels
for(d in c(12,24,36)){
  plot.new()
  text(0.5,0.5,paste(d,"-hour duration", sep=""), font=2, cex=1.2)
}

# Row labels + QQ plots
k <- 1

for(rp in c(5,10,25,50,100)){
  
  # Left-side RP label
  plot.new()
  text(0.5,0.5,paste(rp,"-year RP", sep=""), srt=90, font=2, cex=1.2)
  
  for(dur in c(12,24,36)){
    
    obs_col = paste0("obs_rp", rp)
    sim_col = paste0("sim_rp", rp)
    
    rows = which(df_tc$dur == dur)
    
    q_obs = quantile(df_tc[rows, ][[obs_col]], seq(0.01,0.99,0.01))
    q_sim = quantile(df_tc[rows, ][[sim_col]], seq(0.01,0.99,0.01))
    
    xlim_min = ylim_min = min(q_obs,q_sim)
    xlim_max = ylim_max = max(q_obs,q_sim)
    
    par(mar=c(2,2,2,1))
    
    plot(q_obs, q_sim,
         xlim=c(xlim_min,xlim_max),
         ylim=c(ylim_min,ylim_max),
         xlab="",
         ylab="",
         pch=16,
         cex=0.7)
    
    abline(0,1,col='red',lwd=2)
    
    k <- k + 1
  }
}



# Non_tc
lay <- matrix(c(
  1,  2,  3,  4,
  5,  6,  7,  8,
  9, 10, 11, 12,
  13, 14, 15, 16,
  17, 18, 19, 20,
  21, 22, 23, 24
), nrow = 6, byrow = TRUE)

layout(lay, widths = c(1,4,4,4), heights = c(1,4,4,4,4,4))

par(mar=c(0,0,0,0))

# Empty top-left corner
plot.new()

# Top duration labels
for(d in c(12,24,36)){
  plot.new()
  text(0.5,0.5,paste(d,"-hour duration", sep=""), font=2, cex=1.2)
}

# Row labels + QQ plots
k <- 1

for(rp in c(5,10,25,50,100)){
  
  # Left-side RP label
  plot.new()
  text(0.5,0.5,paste(rp,"-year RP", sep=""), srt=90, font=2, cex=1.2)
  
  for(dur in c(12,24,36)){
    
    obs_col = paste0("obs_rp", rp)
    sim_col = paste0("sim_rp", rp)
    
    rows = which(df_tc$dur == dur)
    
    q_obs = quantile(df_non_tc[rows, ][[obs_col]], seq(0.01,0.99,0.01))
    q_sim = quantile(df_non_tc[rows, ][[sim_col]], seq(0.01,0.99,0.01))
    
    xlim_min = ylim_min = min(q_obs,q_sim)
    xlim_max = ylim_max = max(q_obs,q_sim)
    
    par(mar=c(2,2,2,1))
    
    plot(q_obs, q_sim,
         xlim=c(xlim_min,xlim_max),
         ylim=c(ylim_min,ylim_max),
         xlab="",
         ylab="",
         pch=16,
         cex=0.7)
    
    abline(0,1,col='red',lwd=2)
    
    k <- k + 1
  }
}


# Distribution of factors
hist(df_sim$ntr_factor[df_sim$sample=="ntr_tc"])
length(which(df_sim$ntr_factor[df_sim$sample=="ntr_tc"] >1)) / length(which(df_sim$sample=="ntr_tc")) #90% 
hist(df_sim$ntr_factor[df_sim$sample=="ntr_non_tc"]) 
length(which(df_sim$ntr_factor[df_sim$sample=="ntr_non_tc"] >1)) / length(which(df_sim$sample=="ntr_non_tc")) #76%

# Rad
obs_wl_max_df = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_wl_max_df.csv',sep=""))[,-1]
sim_wl_max_df = read.csv(paste0('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/', name, '_sim_wl_max_df.csv'))[,-1]

par(mfrow=c(3,2))
par(mar=c(4.2,4.5,0.5,0.5))
hist(obs_wl_max_df$wl_max_wl,freq=F)
hist(sim_wl_max_df$wl_max_wl,freq=F)

hist(obs_wl_max_df$wl_max_ntr,freq=F)
hist(sim_wl_max_df$wl_max_ntr,freq=F) 

hist(obs_wl_max_df$wl_max_pred,freq=F) 
hist(sim_wl_max_df$wl_max_pred,freq=F)


#Are highest water levels in ntr sample
# Read IDs corresponding to highest water-level events
wl_ids <- read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name, '_wl_event_data_df.csv', sep = ""))

# Separate NTR peak IDs by storm type
tc_ids <- df_samples$ids_con_ntr[df_samples$type_con_ntr == "tc"]
non_tc_ids <- df_samples$ids_con_ntr[df_samples$type_con_ntr == "non_tc"]

# Count the number of highest water-level events occurring within ±72 h of each NTR peak
tc_index <- sapply(tc_ids, function(x) sum(abs(x - wl_ids$ids) < 72))
non_tc_index <- sapply(non_tc_ids, function(x) sum(abs(x - wl_ids$ids) < 72))

# Proportion of NTR peaks associated with at least one highest water-level event
tc_prop <- mean(tc_index > 0)
non_tc_prop <- mean(non_tc_index > 0)


# Separate NTR peak IDs by storm type
tc_ids <- wl_ids$ids[wl_ids$type == "tc"]
non_tc_ids <- wl_ids$ids[wl_ids$type == "non_tc"]

# Count the number of highest water-level events occurring within ±72 h of each NTR peak
tc_index <- sapply(tc_ids, function(x) sum(abs(x - df_samples$ids_con_ntr) < 72))
non_tc_index <- sapply(non_tc_ids, function(x) sum(abs(x - df_samples$ids_con_ntr) < 72))

# Proportion of NTR peaks associated with at least one highest water-level event
tc_prop <- mean(tc_index > 0)
non_tc_prop <- mean(non_tc_index > 0)



# Read IDs corresponding to highest water-level events
wl_ids <- read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name, '_wl_event_data_df.csv', sep = ""))

# Separate NTR peak IDs by storm type
tc_ids <- df_samples$ids_con_ntr[df_samples$type_con_ntr == "tc"]
non_tc_ids <- df_samples$ids_con_ntr[df_samples$type_con_ntr == "non_tc"]

# Count the number of highest water-level events occurring within ±72 h of each NTR peak
tc_index <- sapply(tc_ids, function(x) sum(abs(x - wl_ids$ids) < 36))
non_tc_index <- sapply(non_tc_ids, function(x) sum(abs(x - wl_ids$ids) < 36))

# Proportion of NTR peaks associated with at least one highest water-level event
tc_prop <- mean(tc_index > 0)
non_tc_prop <- mean(non_tc_index > 0)


# Separate NTR peak IDs by storm type
tc_ids <- wl_ids$ids[wl_ids$type == "tc"]
non_tc_ids <- wl_ids$ids[wl_ids$type == "non_tc"]

# Count the number of highest water-level events occurring within ±72 h of each NTR peak
tc_index <- sapply(tc_ids, function(x) sum(abs(x - df_samples$ids_con_ntr) < 36))
non_tc_index <- sapply(non_tc_ids, function(x) sum(abs(x - df_samples$ids_con_ntr) < 36))

# Proportion of NTR peaks associated with at least one highest water-level event
tc_prop <- mean(tc_index > 0)
non_tc_prop <- mean(non_tc_index > 0)




#Finding minimum differences
#Observed NTR
min(diff(df_samples$ids_con_ntr[order(df_samples$ids_con_ntr)]))
#Observed wls
min(diff(wl_ids$ids[order(wl_ids$ids)]))


png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/wl_events_decomposition.png",width=15,height=7.5,res=300,units='in')
par(mar=c(2,4.5,0.1,0.1))

#Plotting highest wl events - differentiating tide and ntr
event_ids <- wl_ids$ids
n <- length(event_ids)

pred <- cora_df$pred[event_ids]
ntr  <- cora_df$ntr[event_ids]


plot(NA, xaxt='n',
     xlim = c(7, n - 7),
     ylim = c(0, max(pred + ntr, na.rm = TRUE) * 1.1),
     xlab = "",
     ylab = "Water level (m MSL)",
     xaxt = "n", cex.lab=1.5, cex.axis=1.5)

# Shade TC events
for(i in seq_along(event_ids)) {
  if(wl_ids$types[i] == "tc") {
    rect(i - 0.5, par("usr")[3],
         i + 0.5, par("usr")[4],
         col = "grey80",
         border = NA)
  }
}

#
#axis(1, at=1:n, substring(cora_df$date_time[wl_ids$ids],1,8),las=2,cex.axis=0.5)

#Is event within +/- 72 hours of a ntr peak
alpha_val <- ifelse(sapply(wl_ids$ids, function(x) sum(abs(x - df_samples$ids_con_ntr) < 72))>0,1,0.45)
tide_col <- "#4575B4"
ntr_col  <- "#D73027"

# Predicted tide component
rect(seq_along(event_ids) - 0.4,
     0,
     seq_along(event_ids) + 0.4,
     pred,
     col = alpha(tide_col,alpha_val),
     border = NA)

# NTR stacked on top
rect(seq_along(event_ids) - 0.4,
     pred,
     seq_along(event_ids) + 0.4,
     pred + ntr,
     col = alpha(ntr_col,alpha_val),
     border = NA)

#axis(1, at = seq_along(event_ids), labels = event_ids, las = 2, cex.axis = 0.7)
mtext("Event", side=1, line=0.75, cex=1.5)

legend("topright",
       fill = c(tide_col, ntr_col, "grey80"), 
       legend = c("Tide", "NTR", "TC event"), border = NA, cex = 1.2)

dev.off()

#Order by tide
png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/wl_events_decomposition_by_tide.png",width=15,height=7.5,res=300,units='in')
par(mar=c(2,4.5,0.1,0.1))

#Plotting highest wl events - differentiating tide and ntr
ord <- order(cora_df$pred[wl_ids$ids])

event_ids <- wl_ids$ids[ord]
event_types <- wl_ids$types[ord]

n <- length(event_ids)

pred <- cora_df$pred[event_ids]
ntr  <- cora_df$ntr[event_ids]

plot(NA, xaxt='n',
     xlim = c(7, n - 7),
     ylim = c(0, max(pred + ntr, na.rm = TRUE) * 1.1),
     xlab = "",
     ylab = "Water level (m MSL)",
     xaxt = "n", cex.lab=1.5, cex.axis=1.5)

# Shade TC events
for(i in seq_along(event_ids)) {
  if(event_types[i] == "tc") {
    rect(i - 0.5, par("usr")[3],
         i + 0.5, par("usr")[4],
         col = "grey80",
         border = NA)
  }
}

#
#axis(1, at=1:n, substring(cora_df$date_time[wl_ids$ids],1,8),las=2,cex.axis=0.5)

#Is event within +/- 72 hours of a ntr peak
alpha_val <- ifelse(sapply(wl_ids$ids[ord], function(x) sum(abs(x - df_samples$ids_con_ntr) < 72))>0,1,0.45)

tide_col <- "#4575B4"
ntr_col  <- "#D73027"

# Predicted tide component
rect(seq_along(event_ids) - 0.4,
     0,
     seq_along(event_ids) + 0.4,
     pred,
     col = alpha(tide_col,alpha_val),
     border = NA)

# NTR stacked on top
rect(seq_along(event_ids) - 0.4,
     pred,
     seq_along(event_ids) + 0.4,
     pred + ntr,
     col = alpha(ntr_col,alpha_val),
     border = NA)

#axis(1, at = seq_along(event_ids), labels = event_ids, las = 2, cex.axis = 0.7)
mtext("Event", side=1, line=0.75, cex=1.5)

legend("topright",
       fill = c(tide_col, ntr_col, "grey80"), 
       legend = c("Tide", "NTR", "TC event"), border = NA, cex = 1.2)

dev.off()

#Order by ntr
png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/wl_events_decomposition_by_ntr.png",width=15,height=7.5,res=300,units='in')
par(mar=c(2,4.5,0.1,0.1))

#Plotting highest wl events - differentiating tide and ntr
ord <- order(cora_df$ntr[wl_ids$ids])

event_ids <- wl_ids$ids[ord]
event_types <- wl_ids$types[ord]

n <- length(event_ids)

pred <- cora_df$pred[event_ids]
ntr  <- cora_df$ntr[event_ids]

plot(NA, xaxt='n',
     xlim = c(7, n - 7),
     ylim = c(0, max(pred + ntr, na.rm = TRUE) * 1.1),
     xlab = "",
     ylab = "Water level (m MSL)",
     xaxt = "n", cex.lab=1.5, cex.axis=1.5)

# Shade TC events
for(i in seq_along(event_ids)) {
  if(event_types[i] == "tc") {
    rect(i - 0.5, par("usr")[3],
         i + 0.5, par("usr")[4],
         col = "grey80",
         border = NA)
  }
}

#
#axis(1, at=1:n, substring(cora_df$date_time[wl_ids$ids],1,8),las=2,cex.axis=0.5)

#Is event within +/- 72 hours of a ntr peak
alpha_val <- ifelse(sapply(wl_ids$ids[ord], function(x) sum(abs(x - df_samples$ids_con_ntr) < 72))>0,1,0.45)

tide_col <- "#4575B4"
ntr_col  <- "#D73027"

# Predicted tide component
rect(seq_along(event_ids) - 0.4,
     0,
     seq_along(event_ids) + 0.4,
     pred,
     col = alpha(tide_col,alpha_val),
     border = NA)

# NTR stacked on top
rect(seq_along(event_ids) - 0.4,
     pred,
     seq_along(event_ids) + 0.4,
     pred + ntr,
     col = alpha(ntr_col,alpha_val),
     border = NA)

#axis(1, at = seq_along(event_ids), labels = event_ids, las = 2, cex.axis = 0.7)
mtext("Event", side=1, line=0.75, cex=1.5)

legend("topright",
       fill = c(tide_col, ntr_col, "grey80"), 
       legend = c("Tide", "NTR", "TC event"), border = NA, cex = 1.2)

dev.off()


#Simulated events component decomposition
png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/sim_wl_events_decomposition.png",width=15,height=7.5,res=300,units='in')
par(mar=c(2,4.5,0.1,0.1))

tide_col <- "#4575B4"
ntr_col  <- "#D73027"

#Plotting highest wl events - differentiating tide and ntr
event_ids <- max.col(wl$sim_wl_trim) 
n <- length(event_ids)

pred <- wl$sim_pred_trim[cbind(seq_len(nrow(wl$sim_pred_trim)), event_ids)]
ntr  <- wl$sim_ntr_trim[cbind(seq_len(nrow(wl$sim_ntr_trim)), event_ids)]


plot(NA, xaxt='n',
     xlim = c(7, n - 7),
     ylim = c(0, max(pred + ntr, na.rm = TRUE) * 1.1),
     xlab = "",
     ylab = "Water level (m MSL)",
     xaxt = "n", cex.lab=1.5, cex.axis=1.5)

# Shade TC events
for(i in seq_along(event_ids)) {
  if(df_sim$type[i] == "tc") {
    rect(i - 0.5, par("usr")[3],
         i + 0.5, par("usr")[4],
         col = "grey80",
         border = NA)
  }
}

# Predicted tide component
rect(seq_along(event_ids) - 0.4,
     0,
     seq_along(event_ids) + 0.4,
     pred,
     col = tide_col,
     border = NA)

# NTR stacked on top
rect(seq_along(event_ids) - 0.4,
     pred,
     seq_along(event_ids) + 0.4,
     pred + ntr,
     col = ntr_col,
     border = NA)

#axis(1, at = seq_along(event_ids), labels = event_ids, las = 2, cex.axis = 0.7)
mtext("Event", side=1, line=0.75, cex=1.5)

legend("topright",
       fill = c(tide_col, ntr_col, "grey80"),
       legend = c("Predicted", "NTR", "TC event"), border = NA, cex = 1.2)

dev.off()


#Compare GPDs of maximum water levels around ntr_peaks (observed vs simulated)
obs_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_lower.csv',sep=""))[,-1]
obs_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_est.csv',sep=""))[,-1]
obs_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_obs_ntr_wl_peak_72_hr_gpd_combined_upper.csv',sep=""))[,-1]

sim_lower = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_lower.csv',sep=""))[,-1]
sim_est = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_est.csv',sep=""))[,-1]
sim_upper = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_gpd_combined_upper.csv',sep=""))[,-1]

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/gpd_ntr_wl_peak_72_hr_return_level.png",width=10,height=7.5,res=300,units='in')
par(mar=c(4.2,4.5,0.5,0.5))
plot(0, type='n', xlim=c(1,100), ylim=c(1,3.5),
     xaxt='n',yaxt='n',xlab="", ylab="")

mtext("Return Period (years)", side=1,line=2.75)
mtext("WL (m)", side = 2, line = 2.75)

axis(1, at =  c(1,20,40,60,80,100), labels =   c(1,20,40,60,80,100))      
axis(1, at = seq(10,90,20), labels = FALSE, tck= -0.005)  

axis(2, at = seq(1,3.5,0.5), labels = seq(1,3.5,0.5))    
axis(2, at = seq(1.25,3.25,0.25), labels = FALSE, tck=-0.005)  

# ── Obs envelope ──────────────────────────────────────────────
obs_lower_env <- apply(obs_est, 2, quantile, 0.025)
obs_upper_env <- apply(obs_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(obs_upper_env, rev(obs_lower_env)),
        col = adjustcolor("black", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(obs_est, 2, mean), col = "black", lwd = 2)

# ── Sim envelope ──────────────────────────────────────────────
sim_lower_env <- apply(sim_est, 2, quantile, 0.025)
sim_upper_env <- apply(sim_est, 2, quantile, 0.975)

polygon(c(1:100, 100:1),
        c(sim_upper_env, rev(sim_lower_env)),
        col = adjustcolor("red", alpha.f = 0.15),
        border = NA)
lines(1:100, apply(sim_est, 2, mean), col = "red", lwd = 2)

legend("topleft",
       legend = c("Observed", "Simulated"),
       lty    = 1, lwd = 2, col = c("black","red"),
       bty    = "n")


dev.off()

##
obs = sapply(ntr_event_data$ids,function(x){max(cora_df$wl[(x-72):(x+72)])})

# Read IDs corresponding to highest water-level events
wl_ids <- read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name, '_wl_event_data_df.csv', sep = ""))
wl_vals = cora_df$wl[wl_ids$ids]

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/obs_wl_ntr_wl_72_hr.png",width=10,height=7.5,res=300,units='in')

# shared breaks so bins match
brks <- seq(min(c(obs, wl_vals)), max(c(obs, wl_vals)), length.out = 30)

#water level events
hist(obs, breaks = brks, freq = FALSE, col = alpha("black",0.5), border = NA, xlim = range(brks), main = "", xlab = "Water level (m MSL)")

#max water levels within +/- 72 hours of ntr peaks 
hist(wl_vals, breaks = brks, freq = FALSE, col = alpha("red",0.2), border = NA, add = TRUE)

legend("topright", legend = c("NTR events", "WL events"), fill = c(alpha("red",0.5), alpha("black",0.5)), bty = "n")

dev.off()



##Histograms of observed vs simulated at time of water level peaks 

wl_ids <- read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name, '_wl_event_data_df.csv', sep = ""))
wl_vals = cora_df$wl[wl_ids$ids]


png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Histograms_wl_sim_tc.png",width=10,height=7.5,res=300,units='in')
par(mfrow=c(3,2))

#Peaks of tc events
wl_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_tc'),]) 

# Tides
brks <- pretty(range(c(cora_df$pred[wl_ids$ids[wl_ids$types=="tc"]], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)])), n = 20)
hist(cora_df$pred[wl_ids$ids[wl_ids$types=="tc"]],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "Tide (m)")
hist(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "Tide (m)")

# NTR
brks <- pretty(range(c(cora_df$ntr[wl_ids$ids[wl_ids$types=="tc"]], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)])), n = 20)
hist(cora_df$ntr[wl_ids$ids[wl_ids$types=="tc"]],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "NTR (m)")
hist(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)],pch=16,breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "NTR (m)")

#Dependence
plot(pobs(cora_df$pred[wl_ids$ids[wl_ids$types=="tc"]]),pobs(cora_df$ntr[wl_ids$ids[wl_ids$types=="tc"]]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
plot(pobs(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)]),
     pobs(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)]),pch=16,xlab="Tide CDF",ylab="NTR CDF")

dev.off()

png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Histograms_wl_sim_non_tc.png",width=10,height=7.5,res=300,units='in')
par(mfrow=c(3,2))

#peaks of non-tc events
wl_non_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_non_tc'),]) 

#Tide
brks <- pretty(range(c(cora_df$pred[wl_ids$ids[wl_ids$types=="non_tc"]], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)])), n = 20)
hist(cora_df$pred[wl_ids$ids[wl_ids$types=="non_tc"]],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "Tide (m)")
hist(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "Tide (m)")

#NTR
brks <- pretty(range(c(cora_df$ntr[wl_ids$ids[wl_ids$types=="non_tc"]], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)])), n = 20)
hist(cora_df$ntr[wl_ids$ids[wl_ids$types=="non_tc"]],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "NTR (m)")
hist(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "NTR (m)")

#dDependence
plot(pobs(cora_df$pred[wl_ids$ids[wl_ids$types=="non_tc"]]),pobs(cora_df$ntr[wl_ids$ids[wl_ids$types=="non_tc"]]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
plot(pobs(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)]),
     pobs(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
dev.off()



##Histograms of observed vs simulated peak water levels within 72hours of ntr peak 
obs_wl_72_id_tc = sapply(ntr_event_data$ids[ntr_event_data$types=="tc"], function(x) {(x - 72) + which.max(cora_df$wl[(x-72):(x+72)]) - 1})
obs_wl_72_id_non_tc = sapply(ntr_event_data$ids[ntr_event_data$types=="non_tc"], function(x) {(x - 72) + which.max(cora_df$wl[(x-72):(x+72)]) - 1})

#png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Histograms_wl_72_sim_tc.png",width=10,height=7.5,res=300,units='in')
par(mfrow=c(3,2))

#Peaks of tc events
wl_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_tc'),]) 

# Tides
brks <- pretty(range(c(cora_df$pred[obs_wl_72_id_tc], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)])), n = 20)
hist(cora_df$pred[obs_wl_72_id_tc],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "Tide (m)")
hist(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "Tide (m)")

# NTR
brks <- pretty(range(c(cora_df$ntr[obs_wl_72_id_tc], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)])), n = 20)
hist(cora_df$ntr[obs_wl_72_id_tc],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "NTR (m)")
hist(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)],pch=16,breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "NTR (m)")

#Dependence
plot(pobs(cora_df$pred[obs_wl_72_id_tc]),pobs(cora_df$ntr[obs_wl_72_id_tc]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
plot(pobs(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)]),
     pobs(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)]),pch=16,xlab="Tide CDF",ylab="NTR CDF")

#dev.off()

#png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/Histograms_wl_72_sim_non_tc.png",width=10,height=7.5,res=300,units='in')
par(mfrow=c(3,2))

#peaks of non-tc events
wl_non_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_non_tc'),]) 

#Tide
brks <- pretty(range(c(cora_df$pred[obs_wl_72_id_non_tc], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)])), n = 20)
hist(cora_df$pred[obs_wl_72_id_non_tc],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "Tide (m)")
hist(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "Tide (m)")

#NTR
brks <- pretty(range(c(cora_df$ntr[obs_wl_72_id_non_tc], wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)])), n = 20)
hist(cora_df$ntr[obs_wl_72_id_non_tc],breaks = brks,freq=FALSE,ylab="Density", main = "", xlab = "NTR (m)")
hist(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)],breaks = brks,col="red",border = "red",freq=FALSE, ylab="Density", main = "", xlab = "NTR (m)")

#dDependence
plot(pobs(cora_df$pred[obs_wl_72_id_non_tc]),pobs(cora_df$ntr[obs_wl_72_id_non_tc]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
plot(pobs(wl_sim_pred_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)]),
     pobs(wl_sim_ntr_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)]),pch=16,xlab="Tide CDF",ylab="NTR CDF")
#dev.off()






##NTR without fitting GPD by stormtype
# ── Return level comparison: Observed vs Simulated NTR, by storm type ───────

# Function to compute empirical return periods
compute_return_periods <- function(x, n_years) {
  x_sorted <- sort(x, decreasing = TRUE)
  n <- length(x_sorted)
  rank <- 1:n
  rp <- (n_years + 1) / rank   # Weibull plotting position, in years
  data.frame(value = x_sorted, return_period = rp)
}

# Bootstrap CI for return levels
bootstrap_return_levels <- function(sim_data, n_years_sim, n_boot = 100, n,  
                                    rp_grid = c(1, 2, 5, 10, 20, 50, 100)) {
  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(rp_grid))
  
  for (b in 1:n_boot) {
    resample <- sample(sim_data, n, replace = TRUE)
    rp_obj <- compute_return_periods(resample, n_years_sim)
    boot_mat[b, ] <- approx(rp_obj$return_period, rp_obj$value, 
                            xout = rp_grid, rule = 2)$y
  }
  data.frame(
    return_period = rp_grid,
    lower = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE),
    upper = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE),
    median = apply(boot_mat, 2, median, na.rm = TRUE)
  )
}


obs_tc      <- compute_return_periods(cora_df$ntr[ntr_event_data$ids[ntr_event_data$types=="tc"]], n_years = 44)
obs_non_tc  <- compute_return_periods(cora_df$ntr[ntr_event_data$ids[ntr_event_data$types=="non_tc"]], n_years = 44)
sim_tc      <- compute_return_periods(df_sim$ntr[df_sim$sample=="ntr_tc"], n_years = length(df_sim$sample=="tc")/5)
sim_non_tc  <- compute_return_periods(df_sim$ntr[df_sim$sample=="ntr_non_tc"], n_years = length(df_sim$sample=="non_tc")/5)

rp_grid <- seq(0, 46, 0.01)
ci_tc     <- bootstrap_return_levels(df_sim$ntr[df_sim$sample=="ntr_tc"], n_years_sim= 44, rp_grid = rp_grid[(rp_grid > min(obs_tc$return_period)) & (rp_grid < max(obs_tc$return_period))], n=13)
ci_non_tc <- bootstrap_return_levels(df_sim$ntr[df_sim$sample=="ntr_non_tc"], n_years_sim=44, rp_grid = rp_grid[(rp_grid > min(obs_non_tc$return_period)) & (rp_grid < max(obs_non_tc$return_period))], n=44*5-13)


png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/GPD_NTR_WL_3x2.png",width = 6.5, height = 7.5,res=300,units='in')

layout(mat=matrix(c(1:8,9,9),ncol=2,byrow=T),heights=c(0.1,1,1,1,0.2))


par(mar = c(0, 4, 0, 0.25)) 

plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1,1), ylim = c(-1,1))
text(0, 0, "TC", cex = 1.2, font = 2)

# Panel 2: Non-TC header
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1,1), ylim = c(-1,1))
text(0, 0, "Non-TC", cex = 1.2, font = 2)

par(mar=c(4,4,0.25,0.25))
plot(obs_tc$return_period, obs_tc$value, log="x", type = "n",
     xlim = c(min(obs_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", main = "", xaxt='n', yaxt='n')
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak NTR (m)", side=2,line=2.3, cex=0.7)
axis(1, c(5,10,20,50), labels=FALSE)
axis(2, seq(0.5,2.5,0.5), labels=FALSE)
mtext(c(5,10,20,50), at= c(5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(0.5,2.5,0.5), at = seq(0.5,2.5,0.5), side=2, line=0.75, cex=0.7)
polygon(c(ci_tc$return_period, rev(ci_tc$return_period)),
        c(ci_tc$lower, rev(ci_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_tc$return_period, ci_tc$median, col = "red", lwd = 1)
points(obs_tc$return_period, obs_tc$value, pch = 19, col = "black", cex=0.5)
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(a)", cex = 1, adj = c(0,1))


plot(obs_non_tc$return_period, obs_non_tc$value, log = "x", type = "n",
     xlim = c(min(obs_non_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", main = "", xaxt='n', yaxt='n')
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak NTR (m)", side=2,line=2.3, cex=0.7)
axis(1, c(0.2,0.5,1,2,5,10,20,50), labels=FALSE)
axis(2, seq(0.5,2.5,0.5), labels=FALSE)
mtext(c(0.2,0.5,1,2,5,10,20,50), at= c(0.2,0.5,1,2,5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(0.5,2.5,0.5), at = seq(0.5,2.5,0.5), side=2, line=0.75, cex=0.7)
polygon(c(ci_non_tc$return_period, rev(ci_non_tc$return_period)),
        c(ci_non_tc$lower, rev(ci_non_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_non_tc$return_period, ci_non_tc$median, col = "red", lwd = 1)
points(obs_non_tc$return_period, obs_non_tc$value, pch = 19, col = "black",cex=0.5)
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(b)", cex = 1, adj = c(0,1))

#par(mfrow = c(1, 1))
#dev.off()

#Read in simulated water levels
wl_sim_wl_trim = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_wl_trim.csv',sep=""))[,-1]
wl_sim_ntr_trim = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_ntr_trim.csv',sep=""))[,-1]
wl_sim_pred_trim = read.csv(paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/',name,'_sim_wl_pred_trim.csv',sep=""))[,-1]

#Return period plots for max water levels within +/- 72 hours of a ntr peak
obs = sapply(ntr_event_data$ids[ntr_event_data$types=="tc"],function(x){max(cora_df$wl[(x-72):(x+72)])})
obs_tc      <- compute_return_periods(obs, n_years = 44)
obs = sapply(ntr_event_data$ids[ntr_event_data$types=="non_tc"],function(x){max(cora_df$wl[(x-72):(x+72)])})
obs_non_tc  <- compute_return_periods(obs, n_years = 44)

rp_grid <- seq(0, 46, 0.01)
wl_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_tc'),]) 
ci_tc     <- bootstrap_return_levels(wl_sim_wl_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)], n_years_sim= 44, rp_grid = rp_grid[(rp_grid > min(obs_tc$return_period)) & (rp_grid < max(obs_tc$return_period))], n=13)
wl_non_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_non_tc'),]) 
ci_non_tc <- bootstrap_return_levels(wl_sim_wl_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)], n_years_sim=44, rp_grid = rp_grid[(rp_grid > min(obs_non_tc$return_period)) & (rp_grid < max(obs_non_tc$return_period))], n=44*5-13)

#png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/GPD_WL_Peaks_within_72_hours_NTR_Peak.png",width=10,height=7.5,res=300,units='in')

# Add shaded band to panel 1
#par(mfrow = c(1, 2))

plot(obs_tc$return_period, obs_tc$value, log="x", type = "n",
     xlim = c(min(obs_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", xaxt='n', yaxt='n')
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak WL near NTR peak (m)", side=2,line=2.3, cex=0.7)
axis(1, c(5,10,20,50), labels=FALSE)
axis(2, seq(1,3,0.5), labels=FALSE)
mtext(c(5,10,20,50), at= c(5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(1,3,0.5), at = seq(1,3,0.5), side=2, line=0.75, cex=0.7)
polygon(c(ci_tc$return_period, rev(ci_tc$return_period)),
        c(ci_tc$lower, rev(ci_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_tc$return_period, ci_tc$median, col = "red", lwd = 1)
points(obs_tc$return_period, obs_tc$value, pch = 19, col = "black", cex=0.5)
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(c)", cex = 1, adj = c(0,1))
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")

plot(obs_non_tc$return_period, obs_non_tc$value, log = "x", type = "n",
     xlim = c(min(obs_non_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", xaxt='n', yaxt='n')
axis(1, c(0.2,0.5,1,2,5,10,20,50), labels=FALSE)
axis(2, seq(1,3,0.5), labels=FALSE)
mtext(c(0.2,0.5,1,2,5,10,20,50), at= c(0.2,0.5,1,2,5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(1,3,0.5), at = seq(1,3,0.5), side=2, line=0.75, cex=0.7)
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak WL near NTR peak (m)", side=2,line=2.3, cex=0.7)
polygon(c(ci_non_tc$return_period, rev(ci_non_tc$return_period)),
        c(ci_non_tc$lower, rev(ci_non_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_non_tc$return_period, ci_non_tc$median, col = "red", lwd = 1)
points(obs_non_tc$return_period, obs_non_tc$value, pch = 19, col = "black", cex=0.5)
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(d)", cex = 1, adj = c(0,1))
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")



#Return period plots for water levels
obs_tc      <- compute_return_periods(cora_df$wl[wl_event_data$ids[wl_event_data$types=="tc"]], n_years = 44)
obs_non_tc  <- compute_return_periods(cora_df$wl[wl_event_data$ids[wl_event_data$types=="non_tc"]], n_years = 44)

rp_grid <- seq(0, 46, 0.01)
wl_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_tc'),]) 
ci_tc     <- bootstrap_return_levels(wl_sim_wl_trim[cbind(which(df_sim$sample=='ntr_tc'),wl_tc_max_id)], n_years_sim= 44, rp_grid = rp_grid[(rp_grid > min(obs_tc$return_period)) & (rp_grid < max(obs_tc$return_period))], n=13)
wl_non_tc_max_id = max.col(wl_sim_wl_trim[which(df_sim$sample=='ntr_non_tc'),]) 
ci_non_tc <- bootstrap_return_levels(wl_sim_wl_trim[cbind(which(df_sim$sample=='ntr_non_tc'),wl_non_tc_max_id)], n_years_sim=44, rp_grid = rp_grid[(rp_grid > min(obs_non_tc$return_period)) & (rp_grid < max(obs_non_tc$return_period))], n=44*5-13)

#png("C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/results/correlations/GPD_WL_Peaks.png",width=10,height=7.5,res=300,units='in')

# Add shaded band to panel 1
#par(mfrow = c(1, 2))

plot(obs_tc$return_period, obs_tc$value, log="x", type = "n",
     xlim = c(min(obs_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", main = "",xaxt='n', yaxt='n')
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak water level (m MSL)", side=2,line=2.3, cex=0.7)
axis(1, c(1,2,5,10,20,50), labels=FALSE)
axis(2, seq(1,3,0.5), labels=FALSE)
mtext(c(1,2,5,10,20,50), at= c(1,2,5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(1,3,0.5), at = seq(1,3,0.5), side=2, line=0.75, cex=0.7)
polygon(c(ci_tc$return_period, rev(ci_tc$return_period)),
        c(ci_tc$lower, rev(ci_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_tc$return_period, ci_tc$median, col = "red", lwd = 1)
points(obs_tc$return_period, obs_tc$value, pch = 19, col = "black", cex=0.5)
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(e)", cex = 1, adj = c(0,1))
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")

plot(obs_non_tc$return_period, obs_non_tc$value, log = "x", type = "n",
     xlim = c(min(obs_non_tc$return_period), 50), ylim = range(c(obs_tc$value, ci_tc$lower, ci_tc$upper)),
     xlab = "", ylab = "", main = "", xaxt='n', yaxt='n')
mtext("Return period (years)", side=1,line=2.15, cex=0.7)
mtext("Peak water level (m MSL)", side=2,line=2.3, cex=0.7)
axis(1, c(0.2,0.5,1,2,5,10,20,50), labels=FALSE)
axis(2, seq(1,3,0.5), labels=FALSE)
mtext(c(0.2,0.5,1,2,5,10,20,50), at= c(0.2,0.5,1,2,5,10,20,50), side=1, line=0.65, cex=0.7)
mtext(seq(1,3,0.5), at = seq(1,3,0.5), side=2, line=0.75, cex=0.7)
polygon(c(ci_non_tc$return_period, rev(ci_non_tc$return_period)),
        c(ci_non_tc$lower, rev(ci_non_tc$upper)), col = adjustcolor("red", 0.2), border = NA)
lines(ci_non_tc$return_period, ci_non_tc$median, col = "red", lwd = 1)
points(obs_non_tc$return_period, obs_non_tc$value, pch = 19, col = "black", cex=0.5)
usr <- par("usr")
text(x = 10^(usr[1] + 0.02*(usr[2]-usr[1])), y = usr[4] - 0.02*(usr[4]-usr[3]),
     labels = "(f)", cex = 1, adj = c(0,1))
#legend("topleft", legend = c("Observed", "Simulated (median)"),
#       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), bty = "n")

#par(mfrow = c(1, 1))
#dev.off()




#par(mfrow = c(1, 1))
par(mar = c(0, 0, 0.3, 0))
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", legend = c("Observed", "Simulated (median)"),
       col = c("black", "red"), pch = c(19, NA), lty = c(NA, 1), 
       ncol = 2, pt.cex = 0.5, box.lty = 1)

dev.off()
