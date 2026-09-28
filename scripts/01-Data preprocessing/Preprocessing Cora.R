#Load packages
library(MultiHazard)

#watershed name
name= 'Mullica'

#Cora file
cora_path = paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/CONUS/NJ/',name,'/',name,'_cora_centroid_ts.csv',sep="")

#Read in dataframe
cora_wl_df = read.csv(cora_path)[,-1]

#Convert date column to date
start_time <- as.POSIXct("1979-02-01 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2022-12-31 23:00:00", tz = "UTC")
cora_wl_df$date <- seq(start_time, end_time, by = "hour")

#Results vector
cora_detrend_df = cora_wl_df

#Detrend time series using a 30-day moving window
for(i in 1:5){
 cora_detrend_df[,i+1] = Detrend(cora_wl_df[,c(1,i+1)], Method = "window", Window_Width = 24*30, End_Length = 43830, PLOT = TRUE, x_lab = "Date", y_lab = "Water level (m)")
}

#writing results to csv file
write.csv(cora_detrend_df,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_detrend_ts.csv',sep=""))








#Read in  ntrs
cora_path =paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_5_detrend_with_pred.csv',sep="")

#Read in dataframe
cora_wl_df = read.csv(cora_path)[,-1]

#Convert date column to date
start_time <- as.POSIXct("1979-02-01 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2022-12-31 23:00:00", tz = "UTC")
cora_wl_df$date <- seq(start_time, end_time, by = "hour")

#Results vector
wl_decl_df = cora_wl_df


#Delcuster each time series in turn
wl_decl_df$ntr_decl =  Decluster_SW(cora_wl_df[c(4,3)], Window_Width=24*5)$Declustered


write.csv(wl_decl_df,paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/',name,'_cora_centroid_5_detrend_declust_ts.csv',sep=""))

plot(wl_decl_df$date, wl_decl_df$ntr,pch=16)
points(wl_decl_df$date,wl_decl_df$ntr_decl,col=2,pch=16)


#Find tidal segments by month
#How many hours in 18.6 years
yr.18.6 = round(24*18.6*365.25,0)

#File path
cora_path =paste('C:/Users/ro327497/OneDrive - University of Central Florida/Documents/paper/prav/data/cora/',name,'_cora_centroid_3_detrend_with_pred.csv',sep="")

#Read in dataframe
cora_wl_df = read.csv(cora_path)[,-1]

#Convert date column to date
start_time <- as.POSIXct("1979-02-01 00:00:00", tz = "UTC")
end_time <- as.POSIXct("2022-12-31 23:00:00", tz = "UTC")
cora_wl_df$date <- seq(start_time, end_time, by = "hour")

#Results vector
tides_df = cora_wl_df

#Identify tidal peak
local_maximum(cora_wl_df$tide)

