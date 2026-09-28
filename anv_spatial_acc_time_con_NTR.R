#Load packages
library(parallel)
library(ncdf4)
library(MultiHazard)
library(VineCopula)

##Read in data
#watershed name
name= 'Mullica'

#Cora file
cora_path = paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/CORA/',name,'_cora_centroid_detrend_ntr_ts.csv',sep="")
cora_path_decl = paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/CORA/',name,'_cora_centroid_detrend_declust_ntr_ts.csv',sep="")

#AORC subset lat lon
aorc_subset_lat_lon = read.csv(paste('/anvil/projects/x-ees250144/x-rjane/CONUS/NJ/Data/AORC/',name,'_aorc_subset.csv',sep=""))

#years
years = 1979:2024

#Rainfall event duration
rain_dur = 24

# Helper functions (these need to be defined before parallel processing)
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
  tc_indices <- which(type_ids == "tc")
  non_tc_indices <- which(type_ids == "non-tc")
  
  cor_tc <- if(length(tc_indices) > 1) {
    cor.test(precip_data[tc_indices],
             wl_data[tc_indices], method="kendall")
  } else NA
  
  cor_non_tc <- if(length(non_tc_indices) > 1) {
    cor.test(precip_data[non_tc_indices],
             wl_data[non_tc_indices], method="kendall")
  } else NA
  
  return(list("tc_est" = as.numeric(cor_tc$estimate), "tc_pval" =  as.numeric(cor_tc$p.value),
              "non_tc_est" = as.numeric(cor_non_tc$estimate), "non_tc_pval" = as.numeric(cor_non_tc$p.value)))
}


calculate_copula <- function(event_ids, type_ids, precip_data, wl_data) {
  tc_mask <- type_ids == "tc"
  non_tc_mask <- type_ids == "non-tc"

  cop_tc <- if(sum(tc_mask) > 1) {
    tryCatch(
      BiCopSelect(u1 = pobs(precip_data[tc_mask]),
                  u2 = pobs(wl_data[tc_mask])),
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


#Proportions in Zscheischler et al. (2020)
kl_tail_dependence_prop <- function(x1, y1, x2, y2,
                                    u = 0.9, W = 5, risk = c("sum", "min")) {
  risk <- match.arg(risk)

  # Remove NAs
  ok1 <- !is.na(x1) & !is.na(y1)
  ok2 <- !is.na(x2) & !is.na(y2)
  x1 <- x1[ok1]; y1 <- y1[ok1]
  x2 <- x2[ok2]; y2 <- y2[ok2]

  # Transform margins to standard Pareto
  pareto_transform <- function(v) {
    u <- rank(v) / (length(v) + 1)
    p <- 1 / (1 - u)
    return(p)
  }

  X1 <- pareto_transform(x1)
  Y1 <- pareto_transform(y1)
  X2 <- pareto_transform(x2)
  Y2 <- pareto_transform(y2)

  # Compute risk r
  R1 <- if (risk == "sum") X1 + Y1 else pmin(X1, Y1)
  R2 <- if (risk == "sum") X2 + Y2 else pmin(X2, Y2)

  # Determine threshold quantile
  q1 <- quantile(R1, probs = u, na.rm = TRUE)
  q2 <- quantile(R2, probs = u, na.rm = TRUE)

  # Indices of extremes
  idx1 <- which(R1 > q1)
  idx2 <- which(R2 > q2)

  #Binning using angle
  assign_angle_bins <- function(X, Y, W = 5) {
    angle <- atan2(Y, X)  # range [0, pi/2] for positive X, Y
    cuts <- seq(0, pi/2, length.out = W + 1)
    bins <- cut(angle, breaks = cuts, include.lowest = TRUE, labels = FALSE)
    return(bins)
  }

  bins1 <- assign_angle_bins(X1[idx1], Y1[idx1], W)
  bins2 <- assign_angle_bins(X2[idx2], Y2[idx2], W)

  # Compute proportions in bins
  prop1 <- table(bins1) / length(bins1)
  prop2 <- table(bins2) / length(bins2)

  # Fill missing bins with zero
  all_bins <- 1:W
  prop1_full <- prop1
  prop2_full <- prop2
  for (b in all_bins) {
    if (!(b %in% names(prop1_full))) prop1_full[as.character(b)] <- 0
    if (!(b %in% names(prop2_full))) prop2_full[as.character(b)] <- 0
  }

  # Sort proportions by bin
  prop1_vec <- as.numeric(prop1_full[order(as.numeric(names(prop1_full)))])
  prop2_vec <- as.numeric(prop2_full[order(as.numeric(names(prop2_full)))])

  return(list(prop1 = prop1_vec, prop2 = prop2_vec))
}


con_wl_aorc_point_func <- function(i, aorc_subset_lat_lon, years, name, cora_path, cora_path_decl, rain_dur){
  tryCatch({

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
    
    # Decluster events
    aorc_precip_df_decl <- Decluster_S_SW(Data = aorc_precip_df[,1:2],
                                          Window_Width_Sum = rain_dur,
                                          Window_Width = 5*24)
    
    # Process each CORA location
    results_by_cora <- list()

    #Read in ntr data
    cora_all <- read.csv(cora_path)[,-1]
    cora_decl_all <- read.csv(cora_path_decl)[,-1]
    
    for(l in 1:5) {

      #Read water level time series
      cora_ntr_df = cora_all[,c(1,l+1)]
      colnames(cora_ntr_df) <- c("date", "ntr")
      cora_ntr_df$date = as.POSIXct(cora_ntr_df$date)
      
      #Read in declustered water levels
      cora_ntr_decl_df = cora_decl_all[,c(1,l+1)]
      colnames(cora_ntr_decl_df) <- c("date", "ntr")
      cora_ntr_decl_df$date = as.POSIXct(cora_ntr_decl_df$date)
      
      cora_ntr_decl_df <- HURDAT(cora_ntr_decl_df,
                                 lat.loc = aorc_subset_lat_lon$lat[i],
                                 lon.loc = -aorc_subset_lat_lon$lon[i],
                                 rad = 350)
      
      
      # Event identification
      event_counts <- c(5, 4, 3, 2, 1) * 44
      event_data <- lapply(event_counts, function(n) {
        ids <- order(cora_ntr_decl_df$ntr, decreasing = TRUE)[1:n]
        types <- classify_event_types(ids, cora_ntr_decl_df$Name, start_idx = pmax(1,ids - 2*24), end_idx=pmin(nrow(cora_ntr_decl_df), ids + 24))
        list(ids = ids, types = types)
      })

      # Calculate correlations for each event count
      cors_tc <- numeric(5)
      cors_non_tc <- numeric(5)

      cors_pval_tc <- numeric(5)
      cors_pval_non_tc <- numeric(5)
      
      cop_tc_family <- numeric(5)
      cop_non_tc_family <- numeric(5)

      cop_tc_par1 <- numeric(5)
      cop_non_tc_par1 <- numeric(5)

      cop_tc_par2 <- numeric(5)
      cop_non_tc_par2 <- numeric(5)

      kl_tc_p1 <- numeric(5)
      kl_non_tc_p1 <- numeric(5)

      kl_tc_p2 <- numeric(5)
      kl_non_tc_p2 <- numeric(5)

      kl_tc_p3 <- numeric(5)
      kl_non_tc_p3 <- numeric(5)
      
      n_tc = numeric(5)
      n_non_tc = numeric(5)

      for(k in 1:5) {

        evt = event_data[[k]]
        ntr_ext = evt$ids

        v_ntr <- extract_co_max(event_ids=ntr_ext, data=aorc_precip_df_decl$Totals, window=3*24) 

        correlations <- calculate_correlations(evt$ids, evt$types, v_ntr, cora_ntr_df$ntr[ntr_ext])

        cors_tc[k] <- as.numeric(correlations$tc_est)
        cors_non_tc[k] <- as.numeric(correlations$non_tc_est)

        cors_pval_tc[k] <- as.numeric(correlations$tc_pval)
        cors_pval_non_tc[k] <- as.numeric(correlations$non_tc_pval)
        
        copula = calculate_copula(evt$ids, evt$types,v_ntr,cora_ntr_df$ntr[ntr_ext])

        cop_tc_family[k] <- as.numeric(copula$tc[1])
        cop_non_tc_family[k] <- as.numeric(copula$non_tc[1])

        cop_tc_par1[k] <- as.numeric(copula$tc[2])
        cop_non_tc_par1[k] <- as.numeric(copula$non_tc[2])

        cop_tc_par2[k] <- as.numeric(copula$tc[3])
        cop_non_tc_par2[k] <- as.numeric(copula$non_tc[3])

        tc_indices <- evt$types == "tc"
        non_tc_indices <- evt$types == "non-tc"

        n_tc[k] = sum(tc_indices)
        n_non_tc[k] = sum(non_tc_indices)
     
        kl = kl_tail_dependence_prop(x1=v_ntr[tc_indices],
                                     y1=cora_ntr_df$ntr[ntr_ext][tc_indices],
                                     x2=v_ntr[non_tc_indices],
                                     y2=cora_ntr_df$ntr[ntr_ext][non_tc_indices],
                                     u = 0.9, W = 3, risk = "sum")

        kl_tc_p1[k] <- as.numeric(kl$prop1[1])
        kl_non_tc_p1[k] <- as.numeric(kl$prop2[1])

        kl_tc_p2[k] <- as.numeric(kl$prop1[2])
        kl_non_tc_p2[k] <- as.numeric(kl$prop2[2])

        kl_tc_p3[k] <- as.numeric(kl$prop1[3])
        kl_non_tc_p3[k] <- as.numeric(kl$prop2[3])

      }

      # Store results
      res <- data.frame(
        aorc = i, cora = l,
        cor_5_epy_tc = cors_tc[1], cor_4_epy_tc = cors_tc[2], cor_3_epy_tc = cors_tc[3],
        cor_2_epy_tc = cors_tc[4], cor_1_epy_tc = cors_tc[5],
        cor_5_epy_non_tc = cors_non_tc[1], cor_4_epy_non_tc = cors_non_tc[2],
        cor_3_epy_non_tc = cors_non_tc[3], cor_2_epy_non_tc = cors_non_tc[4],
        cor_1_epy_non_tc = cors_non_tc[5],
        cor_pval_5_epy_tc = cors_pval_tc[1], cor_pval_4_epy_tc = cors_pval_tc[2], cor_pval_3_epy_tc = cors_pval_tc[3],
        cor_pval_2_epy_tc = cors_pval_tc[4], cor_pval_1_epy_tc = cors_pval_tc[5],
        cor_pval_5_epy_non_tc = cors_pval_non_tc[1], cor_pval_4_epy_non_tc = cors_pval_non_tc[2],
        cor_pval_3_epy_non_tc = cors_pval_non_tc[3], cor_pval_2_epy_non_tc = cors_pval_non_tc[4],
        cor_pval_1_epy_non_tc = cors_pval_non_tc[5],
        cop_1_epy_tc = cop_tc_family[1], cop_par1_1_epy_tc = cop_tc_par1[1], cop_par2_1_epy_tc = cop_tc_par2[1],
        cop_1_epy_non_tc  = cop_non_tc_family[1],cop_par1_1_epy_non_tc = cop_non_tc_par1[1],cop_par2_1_epy_non_tc = cop_non_tc_par2[1],
        cop_2_epy_tc = cop_tc_family[2], cop_par1_2_epy_tc = cop_tc_par1[2], cop_par2_2_epy_tc = cop_tc_par2[2],
        cop_2_epy_non_tc  = cop_non_tc_family[2],cop_par1_2_epy_non_tc = cop_non_tc_par1[2],cop_par2_2_epy_non_tc = cop_non_tc_par2[2],
        cop_3_epy_tc = cop_tc_family[3], cop_par1_3_epy_tc = cop_tc_par1[3], cop_par2_3_epy_tc = cop_tc_par2[3],
        cop_3_epy_non_tc  = cop_non_tc_family[3],cop_par1_3_epy_non_tc = cop_non_tc_par1[3],cop_par2_3_epy_non_tc = cop_non_tc_par2[3],
        cop_4_epy_tc = cop_tc_family[4], cop_par1_4_epy_tc = cop_tc_par1[4], cop_par2_4_epy_tc = cop_tc_par2[4],
        cop_4_epy_non_tc  = cop_non_tc_family[4],cop_par1_4_epy_non_tc = cop_non_tc_par1[4],cop_par2_4_epy_non_tc = cop_non_tc_par2[4],
        cop_5_epy_tc = cop_tc_family[5], cop_par1_5_epy_tc = cop_tc_par1[5], cop_par2_5_epy_tc = cop_tc_par2[5],
        cop_5_epy_non_tc  = cop_non_tc_family[5],cop_par1_5_epy_non_tc = cop_non_tc_par1[5],cop_par2_5_epy_non_tc = cop_non_tc_par2[5],
        kl_p1_1_epy_tc = kl_tc_p1[1], kl_p2_1_epy_tc = kl_tc_p2[1], kl_p3_1_epy_tc = kl_tc_p3[1],
        kl_p1_1_epy_non_tc= kl_non_tc_p1[1], kl_p2_1_epy_non_tc= kl_non_tc_p2[1], kl_p3_1_epy_non_tc= kl_non_tc_p3[1],
        kl_p1_2_epy_tc = kl_tc_p1[2], kl_p2_2_epy_tc = kl_tc_p2[2], kl_p3_2_epy_tc = kl_tc_p3[2],
        kl_p1_2_epy_non_tc= kl_non_tc_p1[2], kl_p2_2_epy_non_tc= kl_non_tc_p2[2], kl_p3_2_epy_non_tc= kl_non_tc_p3[2],
        kl_p1_3_epy_tc = kl_tc_p1[3], kl_p2_3_epy_tc = kl_tc_p2[3], kl_p3_3_epy_tc = kl_tc_p3[3],
        kl_p1_3_epy_non_tc= kl_non_tc_p1[3], kl_p2_3_epy_non_tc= kl_non_tc_p2[3], kl_p3_3_epy_non_tc= kl_non_tc_p3[3],
        kl_p1_4_epy_tc = kl_tc_p1[4], kl_p2_4_epy_tc= kl_tc_p2[4], kl_p3_4_epy_tc = kl_tc_p3[4],
        kl_p1_4_epy_non_tc= kl_non_tc_p1[4], kl_p2_4_epy_non_tc= kl_non_tc_p2[4], kl_p3_4_epy_non_tc= kl_non_tc_p3[4],
        kl_p1_5_epy_tc= kl_tc_p1[5], kl_p2_5_epy_tc = kl_tc_p2[5], kl_p3_5_epy_tc = kl_tc_p3[5],
        kl_p1_5_epy_non_tc = kl_non_tc_p1[5], kl_p2_5_epy_non_tc = kl_non_tc_p2[5], kl_p3_5_epy_non_tc = kl_non_tc_p3[5],
        n_5_epy_tc = n_tc[1], n_5_epy_non_tc  = n_non_tc[1],
        n_4_epy_tc = n_tc[2], n_4_epy_non_tc  = n_non_tc[2],
        n_3_epy_tc = n_tc[3], n_3_epy_non_tc  = n_non_tc[3],
        n_2_epy_tc = n_tc[4], n_2_epy_non_tc  = n_non_tc[4],
        n_1_epy_tc = n_tc[5], n_1_epy_non_tc  = n_non_tc[5]
      )

      results_by_cora[[l]] <- res
    }

    return(results_by_cora)

  }, error = function(e) {
    message(paste("Error processing location", i, ":", e$message))
    traceback()  # Add this
    return(NULL)
  })
}

# Bkeene
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
clusterExport(cl, c("cora_path", "cora_path_decl", "aorc_subset_lat_lon", "years", "name",
                    "con_wl_aorc_point_func", "classify_event_types",
                    "calculate_correlations","calculate_copula",
                    "kl_tail_dependence_prop",
                    "rain_dur", "extract_co_max"))

cat("Cluster created")

# Load required packages on workers
clusterEvalQ(cl, {
  library(ncdf4)
  library(MultiHazard)
  library(VineCopula)
})

cat("Loaded required packages on workers")

# Run parallel processing
all_results <- parLapply(cl, 2001:nrow(aorc_subset_lat_lon), con_wl_aorc_point_func,
                         aorc_subset_lat_lon = aorc_subset_lat_lon,
                         cora_path = cora_path, cora_path_decl = cora_path_decl,
                         name=name, rain_dur=rain_dur,
                         years = years)

# Cleanup
stopCluster(cl)

print(all_results)

# Flatten the nested list structure
flattened_list <- unlist(all_results, recursive = FALSE)

# Combine all dataframes into one
combined_df <- do.call(rbind, flattened_list)

# Reset row names to be sequential
rownames(combined_df) <- NULL

# Write to CSV
write.csv(combined_df, paste('/anvil/projects/x-ees250144/x-rjane/CONUS/', name,'_anv_spatial_acc_time_con_ntr_res_detrend_2500.csv', sep=""), row.names = FALSE)
