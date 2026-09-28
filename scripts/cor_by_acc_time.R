
#Load packages
library(parallel)

#watershed name
name= 'Mullica'

#Cora file
cora_path = paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/CORA/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")
cora_path_decl = paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/CORA/',name,'_cora_centroid_3_detrend_declust_ts.csv',sep="")

#AORC subset lat lon
aorc_subset_lat_lon = read.csv(paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/AORC/',name,'_aorc_subset.csv',sep=""))

#years
years = 1979:2024

#Rainfall accumulation times
rain_dur_seq = seq(1,48,1)

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



acc_cor <- function(i, aorc_subset_lat_lon, years, cora_path, cora_path_decl, name, rain_dur_seq){

  tryCatch({
    
#Results vectors
cor_r = numeric(length(rain_dur_seq)) 
cor_r_tc = numeric(length(rain_dur_seq))
cor_r_non_tc = numeric(length(rain_dur_seq))
cor_ntr = numeric(length(rain_dur_seq))
cor_ntr_tc = numeric(length(rain_dur_seq))
cor_ntr_non_tc = numeric(length(rain_dur_seq))   
cor_pval_r = numeric(length(rain_dur_seq))
cor_pval_r_tc = numeric(length(rain_dur_seq))
cor_pval_r_non_tc = numeric(length(rain_dur_seq))
cor_pval_ntr = numeric(length(rain_dur_seq))  
cor_pval_ntr_tc = numeric(length(rain_dur_seq))  
cor_pval_ntr_non_tc = numeric(length(rain_dur_seq)) 
n_r_tc = numeric(length(rain_dur_seq))
n_r_non_tc = numeric(length(rain_dur_seq)) 
n_ntr_tc = numeric(length(rain_dur_seq))  
n_ntr_non_tc  = numeric(length(rain_dur_seq)) 

# AORC data processing
aorc_precip <- numeric()
aorc_time <- list()

for(j in seq_along(years)) {
  
  aorc_file <- nc_open(paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/AORC/',name,'_aorc_sub_', years[j], '.nc',sep=""))
  
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

#Decluster rainfall volume
aorc_precip_df <- HURDAT(aorc_precip_df,
                         lat.loc = aorc_subset_lat_lon$lat[i],
                         lon.loc = -aorc_subset_lat_lon$lon[i],
                         rad = 350)

#Read water level time series
cora_ntr_df = read.csv(cora_path)[,c(1,4)]
cora_ntr_df[,1] =  seq(as.POSIXct("1979-02-01 00:00:00", tz = "UTC"),as.POSIXct("2022-12-31 23:00:00", tz = "UTC"),by = "hour")
colnames(cora_ntr_df) <- c("date", "ntr")
#cora_ntr_df$date = as.POSIXct(cora_ntr_df$date)

#Read in declustered water levels
cora_ntr_decl_df = read.csv(cora_path_decl)[,-1][,c(4,5)]
colnames(cora_ntr_decl_df) <- c("date", "ntr")
cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)

cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                           lat.loc = aorc_subset_lat_lon$lat[i],
                           lon.loc = -aorc_subset_lat_lon$lon[i],
                           rad = 350)

for(k in 1:length(rain_dur_seq)){


  # Decluster events
  aorc_precip_df_decl <- Decluster_S_SW(Data = aorc_precip_df[,1:2],
                                        Window_Width_Sum = rain_dur_seq[k],
                                        Window_Width = 5*24)
  
  # Event identification
  ids <- order(aorc_precip_df_decl$Declustered, decreasing = TRUE)[1:(5*44)]
  types <- classify_event_types(ids, aorc_precip_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(aorc_precip_df), ids + 24))
  rain_event_data = list(ids = ids, types = types)
  
  #Event identification for ntr
  
  ids <- order(cora_ntr_decl_df$ntr, decreasing = TRUE)[1:(5*44)]
  types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 24))
  ntr_event_data = list(ids = ids, types = types)
  
  evt = ntr_event_data
  ntr_ext = evt$ids
  
  #Generate sample con. on ntr
  v_ntr <- extract_co_max(ntr_ext, aorc_precip_df_decl$Totals, window = 3*24)
  
  #Correlation without stratification by generating mechanism
  cor = cor.test(v_ntr,cora_ntr_df$ntr[ntr_ext], method="kendall")
  cor_ntr[k] <- cor$estimate
  cor_pval_ntr[k] <- cor$p.value
  
  #Characteristics of sample con. on ntr
  correlations <- calculate_correlations_with_sig(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])
  print(correlations)
  cor_ntr_tc[k] <- as.numeric(correlations$tc)
  cor_ntr_non_tc[k] <- as.numeric(correlations$non_tc)
  
  cor_pval_ntr_tc[k] <- as.numeric(correlations$pval_tc)
  cor_pval_ntr_non_tc[k] <- as.numeric(correlations$pval_non_tc)
  
  tc_indices <- evt$types == "tc"
  non_tc_indices <- evt$types == "non-tc"
  
  n_ntr_tc[k] = sum(tc_indices)
  n_ntr_non_tc[k] = sum(non_tc_indices)
  
  
  #Generate sample con. on r
  evt <- rain_event_data
  ntr_vals <- extract_co_max(evt$ids, cora_ntr_df$ntr, window = 3*24)
  
  #Correlation without stratification by generating mechanism
  cor = cor.test(aorc_precip_df_decl$Declustered[evt$ids],ntr_vals, method="kendall")
  cor_r[k] <- cor$estimate
  cor_pval_r[k] <- cor$p.value
  
  #Characteristics of sample con. on r
  correlations <- calculate_correlations_with_sig(evt$ids, evt$types,
                                                  aorc_precip_df_decl$Declustered[evt$ids],
                                                  ntr_vals)
  
  cor_r_tc[k] <- as.numeric(correlations$tc)
  cor_r_non_tc[k] <- as.numeric(correlations$non_tc)
  
  cor_pval_r_tc[k] <- as.numeric(correlations$pval_tc)
  cor_pval_r_non_tc[k] <- as.numeric(correlations$pval_non_tc)
  
  tc_r_indices <- evt$types == "tc"
  non_tc_r_indices <- evt$types == "non-tc"
  
  n_r_tc[k] = sum(tc_r_indices)
  n_r_non_tc[k] = sum(non_tc_r_indices)
  
 }

 # Store results
 res <- data.frame(aorc = rep(i,length(rain_dur_seq)),rain_dur_seq, cora = rep(3,length(rain_dur_seq)),
                   cor_r = cor_r, cor_ntr = cor_ntr,
                   cor_r_tc = cor_r_tc, cor_r_non_tc = cor_r_non_tc,
                   cor_ntr_tc = cor_ntr_tc, cor_ntr_non_tc = cor_ntr_non_tc,
                   cor_pval_r = cor_pval_r, cor_pval_ntr = cor_pval_ntr,
                   cor_pval_r_tc = cor_pval_r_tc, cor_pval_r_non_tc = cor_pval_r_non_tc,
                   cor_pval_ntr_tc = cor_pval_ntr_tc, cor_pval_ntr_non_tc = cor_pval_ntr_non_tc,
                   n_r_tc = n_r_tc, n_r_non_tc = n_r_non_tc,
                   n_ntr_tc = n_ntr_tc, n_ntr_non_tc  = n_ntr_non_tc)
 
 return(res)

  }, error = function(e) {
    message(paste("Error processing location", i, ":", e$message))
    traceback()  # Add this
    return(NULL)
  })
}

cpus_per_task <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
cat("Creating cluster with", cpus_per_task, "CPUS\n")
if(is.na(cpus_per_task) || cpus_per_task < 1) {
  cpus_per_task <- 1
  cat("Warning: defaulting to 1 CPU")
} else {
  cat("SLURM_CPUS_PER_TASK =", cpus_per_task)
}
cl <- makeCluster(cpus_per_task, type = "FORK")

# Export all necessary objects
clusterExport(cl, c("cora_path", "cora_path_decl", "aorc_subset_lat_lon", "years", "name","rain_dur_seq",
                    "extract_co_max",
                    "classify_event_types", 
                    "calculate_correlations", "calculate_correlations_with_sig"))

# Load required packages on workers
clusterEvalQ(cl, {
  library(ncdf4)
  library(MultiHazard)
  library(VineCopula)
  library(MASS)
  library(tweedie)
  library(gamlss)
  library(gamlss.mx)
  library(VGAM)
  library(stats4)
  library(copula)
  library(truncnorm)
})

# Run parallel processing
all_results <- parLapply(cl, 2301:2355, acc_cor,
                         aorc_subset_lat_lon = aorc_subset_lat_lon,
                         cora_path = cora_path,
                         cora_path_decl = cora_path_decl,
                         name= name,
                         years = years,
                         rain_dur_seq = rain_dur_seq)

# Cleanup
stopCluster(cl)

print(all_results)

# Flatten the nested list structure
#flattened_list <- unlist(all_results, recursive = FALSE)

# Combine all dataframes into one
combined_df <- do.call(rbind, all_results)

# Reset row names to be sequential
rownames(combined_df) <- NULL

# Write to CSV
write.csv(combined_df, paste('/anvil/projects/x-ees250144/x-rjane/CONUS/',name,'_cor_by_acc_time_2350.csv',sep=""), row.names = FALSE)
