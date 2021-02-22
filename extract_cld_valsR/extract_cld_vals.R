library(date)
library(reshape2)
library(data.table)
library(dplyr)
library(rgdal)
library(raster)

# all strings of letters coming in will remain character
options(stringsAsFactors=FALSE)

master <- read.csv("D:/Landsat_ESPA/NPWRC_DATA_16FEB2018.csv", sep = ",")

#What year? Please change this to the year you want to run
y <- '2003'

#path row as int
master$PATH <- as.integer(master$PATH)
master$DATE <- as.Date(master$DATE)

pathrow <- read.csv("D:/Landsat_ESPA/path row x wetland x chamber x year.csv", sep = ",")

#get the number of unique vals for the year we are running.
pathrow2 <- pathrow[ which(pathrow$YEAR.x == y), ]

#number of path row combos:
pr_combos <- unique(pathrow2$PATH_ROW)
pr_num <- as.numeric(length(unique(pathrow2$PATH_ROW)))

# Setwd to scenes dir
setwd("D:/Landsat_ESPA/2003")

#Pixel QA Files
pixelQA.files <- list.files(pattern = '\\pixel_qa.tif$')

satellite <- substr(pixelQA.files,1,4)
processing <- substr(pixelQA.files,6,9)
PATH <- substr(pixelQA.files,12,13)
ROW <- substr(pixelQA.files,15,16)
year <- substr(pixelQA.files,18,21)
mmonth <- substr(pixelQA.files,22,23) #extra m stays in there to avoid variable call error within r
day <- substr(pixelQA.files,24,25)
cx <- substr(pixelQA.files,36,37)
tx <- substr(pixelQA.files,39,40)
product <- substr(pixelQA.files,42,46)
band <- substr(pixelQA.files,48,49)
extension <- substr(pixelQA.files,51,53)
file <- pixelQA.files

df.qa <- melt(data.frame(satellite, processing, PATH, ROW, year, mmonth,
                         day, cx, tx, product, band, extension, file, stringsAsFactors = F))

#es1df <- melt(data.frame(file, PATH, ROW, year, mmonth, day, stringsAsFactors = F)) #essentials for matching against master dataframe
es1df <- melt(data.frame(file,satellite,PATH,ROW, year, mmonth, day, stringsAsFactors = F)) #essentials for matching against master dataframe
es1df$PATH <- paste(es1df$PATH,es1df$ROW,sep="_")
#
#make a correct combined date column in data
es1df$r_date <- paste(es1df$year,es1df$mmonth,es1df$day, sep = "-")

#What year we looking at? Subset for year
master2 <- master[ which(master$YEAR == y), ]


######BEFORE DATE

es1df$r_date <- as.Date(es1df$r_date) #make a date
es1df <- es1df[order(es1df$r_date),] #reorder

file<-"NULL"
master2$days_before<-NA
master2$ndvi_before<-NA
#z <- as.Date(c(es1df$r_date)) # make it a vector

##i loop is for observations 
##j loop is for pathrows
##k loop is for date differences

for(i in 1:nrow(master2)){
  
  prmas<-pathrow2[which(master2$WETLAND_ID[i]==pathrow2$WETLAND_ID & master2$CHAMBER_ID[i]==pathrow2$chamber),4]
  
  ndvi_mat<-matrix(NA,nrow=length(prmas))
  date_mat<-matrix(NA,nrow=length(prmas))
  
  for(j in 1:length(prmas)){
    
    fly_dates<-es1df$r_date[es1df$PATH==prmas[j]]
    
    if(master2$DATE[i]<min(fly_dates)){break}
    
    dates_bef<- rev(fly_dates[fly_dates<=master2$DATE[i]])
    
    for(k in 1:length(dates_bef)){
      
      file<-es1df[es1df$PATH==prmas[j] & es1df$r_date==dates_bef[k],1]
      b <- brick(file) #store the cloud raster (pixel.qa) as a brick file format
      qa<-extract(b, SpatialPoints(cbind(master2$LONG_CH[i], master2$LAT_CH[i]),proj4string = CRS("+proj=longlat +datum=WGS84"))) #extract the pixel value based on the chamber lat long using the WGS84 datum
      
      if((qa==66 | qa==68| qa==132| qa==130|qa==322 | qa==386| qa==834| qa==898|qa==1346 | qa==324| qa==388| qa==836| qa==900| qa==1348) & !is.na(qa)){break}
      
    }
    
    if(is.na(qa)){
      ndvi_mat[j,]<-NA
      date_mat[j,]<-NA  
    }else{
    
      if((qa==66 | qa==68| qa==132| qa==130|qa==322 | qa==386| qa==834| qa==898|qa==1346 | qa==324| qa==388| qa==836| qa==900| qa==1348)){
      
        ndvi_file<- sub('pixel_qa.tif', '', file) #remove tag
        ndvi_file <- paste(ndvi_file, 'sr_ndvi.tif', sep = '') #add in ndvi tag
    
        b <- brick(ndvi_file) #store the sr ndvi as a brick file format
        ndvi <- extract(b, SpatialPoints(cbind(master2$LONG_CH[i], master2$LAT_CH[i]),proj4string = CRS("+proj=longlat +datum=WGS84"))) 
    
        ndvi_mat[j,]<-ndvi
        date_mat[j,]<-as.numeric(dates_bef[k]-master2$DATE[i])
      }else{
        ndvi_mat[j,]<-NA
        date_mat[j,]<-NA}}
    
    }
  
  if(all(is.na(ndvi_mat))){
    master2$ndvi_before[i]<-NA
    master2$days_before[i]<-NA
    
  }else{
    master2$ndvi_before[i]<-mean(ndvi_mat[which(abs(date_mat)==min(abs(date_mat),na.rm=TRUE)),],na.rm=TRUE)*0.0001
    master2$days_before[i]<-mean(date_mat[which(abs(date_mat)==min(abs(date_mat),na.rm=TRUE)),],na.rm=TRUE)
  }
  
}



####AFTER
file<-"NULL"
master2$days_after<-NA
master2$ndvi_after<-NA
#z <- as.Date(c(es1df$r_date)) # make it a vector

##i loop is for observations 
##j loop is for pathrows
##k loop is for date differences

for(i in 1:nrow(master2)){
  
  prmas<-pathrow2[which(master2$WETLAND_ID[i]==pathrow2$WETLAND_ID & master2$CHAMBER_ID[i]==pathrow2$chamber),4]
  
  ndvi_mat<-matrix(NA,nrow=length(prmas))
  date_mat<-matrix(NA,nrow=length(prmas))
  
  for(j in 1:length(prmas)){
    
    fly_dates<-es1df$r_date[es1df$PATH==prmas[j]]
    
    if(master2$DATE[i]>max(fly_dates)){break}
    
    dates_aft<- (fly_dates[fly_dates>=master2$DATE[i]])
    
    for(k in 1:length(dates_aft)){
      
      file<-es1df[es1df$PATH==prmas[j] & es1df$r_date==dates_aft[k],1]
      b <- brick(file) #store the cloud raster (pixel.qa) as a brick file format
      qa<-extract(b, SpatialPoints(cbind(master2$LONG_CH[i], master2$LAT_CH[i]),proj4string = CRS("+proj=longlat +datum=WGS84"))) #extract the pixel value based on the chamber lat long using the WGS84 datum
      
      if((qa==66 | qa==68| qa==132| qa==130|qa==322 | qa==386| qa==834| qa==898|qa==1346 | qa==324| qa==388| qa==836| qa==900| qa==1348) & !is.na(qa)){break}
      
    }
    
    if(is.na(qa)){
      ndvi_mat[j,]<-NA
      date_mat[j,]<-NA  
    }else{
      
      if((qa==66 | qa==68| qa==132| qa==130|qa==322 | qa==386| qa==834| qa==898|qa==1346 | qa==324| qa==388| qa==836| qa==900| qa==1348)){
        
        ndvi_file<- sub('pixel_qa.tif', '', file) #remove tag
        ndvi_file <- paste(ndvi_file, 'sr_ndvi.tif', sep = '') #add in ndvi tag
        
        b <- brick(ndvi_file) #store the sr ndvi as a brick file format
        ndvi <- extract(b, SpatialPoints(cbind(master2$LONG_CH[i], master2$LAT_CH[i]),proj4string = CRS("+proj=longlat +datum=WGS84"))) 
        
        ndvi_mat[j,]<-ndvi
        date_mat[j,]<-as.numeric(dates_aft[k]-master2$DATE[i])
      }else{
        ndvi_mat[j,]<-NA
        date_mat[j,]<-NA}}
    
  }
  
  if(all(is.na(ndvi_mat))){
    master2$ndvi_after[i]<-NA
    master2$days_after[i]<-NA
    
  }else{
    master2$ndvi_after[i]<-mean(ndvi_mat[which(abs(date_mat)==min(abs(date_mat),na.rm=TRUE)),],na.rm=TRUE)*0.0001
    master2$days_after[i]<-mean(date_mat[which(abs(date_mat)==min(abs(date_mat),na.rm=TRUE)),],na.rm=TRUE)
  }
  
}

master2$ndvi_interpol=master2$ndvi_after-(((master2$ndvi_after-master2$ndvi_before)/(master2$days_after-master2$days_before))*master2$days_after)

write.csv(master2,"D:/Landsat_ESPA/2005/master2003.csv",row.names=F, na="")

