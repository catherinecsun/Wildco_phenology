##############################################################
# Timelapse_Phenopix_3_PhenoFit_Annually.R
# combine VI's across deployment periods/seasons for a site to fit curves to do full 1 year of phenology extraction  
# this does 3 two-stagefits: Spline+TRS, Spline+Klosterman, and Klosterman+Klosterman
# The first stage fits a curve to reduce the influence of a single obs to an cature the seasonal behavir
# the second stage extracts the phenology dates. 
# BE VERY SKEPTICAL of results if the timelapse images do not capture the majority of at least one single growing season
# because the curves are fit to the given data; it does not infer missing data
# things are relative to the peak, 
# soif you dont have data around the peaks, then you cant trust the estimates.
# if you have 'enough' of a growing season, you can trust some of the dates.
##############################################################

####--- load packages ####
library(phenopix)
library(jpeg)
library(rasterImage)
library(lubridate)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(zoo)
library(data.table)

###---- Define and Set things  THIS NEEDS TO BE SPECIFIC TO YOUR PROJECT AND YOUR FOLDER STRUCTURE. ####
OutputDir <- "Output"
setwd(OutputDir)

# if a year at a site has less data than this, the fitted curves and therefore pheno dates are probably crap
# threshold below which you dont want the results for a site in a year
# for the purposes of this demo, this threshold is very low so that we can complete this script with the small batch of timelapse photos
# if you have too few timelapse images, DO NOT TRUST THE CURVES OR EXTRACTED DATES.
minPropOfYear<-0.05 

#the site(s) you want to do, can be a single site or vector of multiple.
Sites<-c("CATH08")

#if your data do not capture whole growing seasons
# how do you want to deal with that?
# No means just still go by the dates of the timelapse and separate the dates according to the real Jan 1 - Dec 31 dates
# Yes means you want to consider a year's worth of data at a time. so that the first image is temporarily considered Jan 1 and then changed back
# you can use Yes (but dont have to) as long as each sequential year's worth of data still captures a whole growing season
DataShift<-"No"#/Yes or No

###--- create a custom trs function ####
# adapted from Cam McLelland's (Coops lab)
cam_trs_cat <- function(spline_data,  percent){ 
  
  # spline_data = seqence of EVI or GCC values 
  # last_doy = total number of values in spline sequence
  # Percent = what percentile you wish to extract 0 to 1
  # doy_seq = julian date associated with each EVI or GCC value
  last_doy<-length(spline_data)
  
  maxval_loc1 <- which(spline_data == max(spline_data, na.rm = TRUE)) # finds the max value location
  maxval_loc <- maxval_loc1[1] # takes the first max value only in the case on of more than 1
  max_val <- max(spline_data, na.rm = TRUE) #finds the max value
  
  first_min_loc1 <- which(spline_data == min(spline_data[1:maxval_loc], na.rm = TRUE)) # finds min before max
  first_min_loc <- first_min_loc1[1]
  first_min <- spline_data[first_min_loc] # before max, min value
  
  last_min_loc1 <- which(spline_data == min(spline_data[maxval_loc:last_doy], na.rm = TRUE)) # finds min after max
  last_min_loc <- last_min_loc1[1]
  last_min <- spline_data[last_min_loc] # after max, min value
  
  sos_val <- as.numeric((max_val - first_min) * percent + first_min )# calculates sos percent value
  eos_val <- as.numeric((max_val - last_min) * percent + last_min )# calculates eos percent value
  
  closest<-function(xv,sv){
    xv[which(abs(xv-sv)==min(abs(xv-sv), na.rm = TRUE))] } # function for finding the closest value associated with eos_val and sos_val
  # the closets function assumes a monotonically increasing/decreasing function, which probably works for the steep sides
  sos_bva <- as.numeric(closest(spline_data[first_min_loc:maxval_loc], sos_val)) # finds the closest value to sos_val
  sos_bva <- sos_bva[1]
  sos_loc <- which(spline_data[first_min_loc:maxval_loc] == sos_bva) # determines where sos_bva is located in the spline data
  sos_loc <- spline_data[first_min_loc:maxval_loc][sos_loc]# cam previously had it as sos_loc[1]
  
  eos_bva <- as.numeric(closest(spline_data[maxval_loc:last_min_loc], eos_val))# finds the closest value to eos_val
  eos_bva <- eos_bva[1]
  eos_loc <- which(spline_data[maxval_loc:last_min_loc] == eos_bva) # determines where eos_bva is located in the spline data
  eos_loc <- spline_data[maxval_loc:last_min_loc][eos_loc]# cam previously had it as eos_loc[1]
  #eos_loc <- maxval_loc + eos_loc # add max value location eos location, not necessary when zoo
  
  #sos <- doy_seq[sos_loc] #determines the doy that is associated with sos_loc
  #eos <- doy_seq[eos_loc] #determines the doy that is associated with eos_loc
  #pop <- doy_seq[maxval_loc] #determines the doy of the max value POP = seasonal max
  out <- data.frame(Event=c(paste0("SoS_",percent),"PoS",paste0("EoS_",percent)),FakeDate=c(index(sos_loc),index(spline_data[maxval_loc1]),index(eos_loc)),
                    Value=c(as.numeric(sos_loc),as.numeric(spline_data[maxval_loc1]),as.numeric(eos_loc)))
  rownames(out)<-NULL
  return(out)
}


###---Do the Fitting/Analysis  #####
for(i in 1:length(Sites)){
  print(Sites[i])
  #per site
  #create a site-specific folder to put the site specific VI results in
  dir.create(paste0("./FitAndPheno_plots/",Sites[i]))
  
  #### get the csv files of VI in the VI folder and read them in; 
  #### 1 list item per deployment period/season 
  VIs_forSite <- intersect(list.files(paste0("./VI/",Sites[i]),pattern=Sites[i]),
                           list.files(paste0("./VI/",Sites[i],"/"),pattern=".csv"))
  VIs_forSite <- unlist(lapply(paste0("./VI/",Sites[i],"/"),paste0,VIs_forSite))
  VIs_forSite <- lapply(VIs_forSite, read.csv)
  
  
  #### prep a version of VI data where the list items are by ROI rather than season 
  nroi<-as.numeric(length(unique(VIs_forSite[[1]]$ROI)))
  VIs_forSite_perROI<- sapply(sapply("ROI",paste0,seq(1:nroi)),function(x) NULL)
  
  ####do the analysis/
  for(j in 1:nroi){ #per roi for that site
    print(paste0("ROI: ",j))
    #subset list_perSeaon to that ROI
    tmp<-lapply(VIs_forSite, subset, ROI==j)
    tmp<-do.call("rbind",tmp)
    VIs_forSite_perROI[[j]] <-tmp
    VIs_forSite_perROI[[j]]$date <-as.POSIXct(VIs_forSite_perROI[[j]]$date)
    DateCol<-which(colnames( VIs_forSite_perROI[[j]])=="date")
    VIs_forSite_perROI[[j]]<-arrange(VIs_forSite_perROI[[j]],date)
    
    #create a lookup index table for date and day within sequence 
    # greenExplore and greenProcess currently use hydro arg to use either Jan1 or Oct1 as the first DOY, cutting out data. STUPID AS.
    # so table look up also creates 'fake dates'
    # this preiviously accomoodated attempts with Algar data to shift the dates of timelapse images
    # because phenopix would otherwise force phenophases into the few months that the end of the year that started the dataset 
    # so that the first image was for Jan 1 to prevent  greenProcess() from subsetting  data if we wanted to a year's *worth* of data
    # ( by hacking it so that the DOY/Index is based on the fake values in the lookup table,so all data points are used)
    # *** THIS FAKE DATE APROACH IS ONLY OK AND SHOULD BE USED ONLY  IF YOU STILL CAPTURE THE ENTIRE GROWING SEASON
    
    min(VIs_forSite_perROI[[j]]$date)
    max(VIs_forSite_perROI[[j]]$date)
    FakeStart<-as.POSIXct(paste0(year(min(VIs_forSite_perROI[[j]]$date)),"-01-01"),format="%Y-%m-%d")
    DateLookup<-data.frame(DateReal=seq.Date(from=as.Date(min(VIs_forSite_perROI[[j]]$date)),
                                             to=as.Date(max(VIs_forSite_perROI[[j]]$date)),
                                             by="day"))
    DateLookup$DayOfSequence <- seq(1,nrow(DateLookup),by=1)
    DateLookup$doy <- as.numeric(format(DateLookup$DateReal, "%j"))
    DateLookup$DateFake<- as.Date(seq(from=FakeStart,by="day",length.out = nrow(DateLookup)))
    DateLookup$doyFake<-as.numeric(format(DateLookup$DateFake, "%j"))
    #DateLookup
    
    #filter and fit the datalength(Years)
    # per list item in VIs_forSite_perROI (ie the VI data per ROI):
    dn <-c('ri.av','gi.av','bi.av') #cols with RGB digital numbers
    brt<- "bri.av" # col with brightness
    
    Filtered.Data <- autoFilter(as.data.frame(VIs_forSite_perROI[[j]]),raw.dn = FALSE, na.fill = TRUE, 
                                filter = c("night", "spline", "max"),
                                filter.options = NULL, plot= TRUE,
                                dn=dn,
                                brt = brt)
    
    # tmp is identical to DayLookup, unless there were days that malfunctioned (there will be fewer rows, but ther indices will be correct)
    tmp_DayLookup<-merge(DateLookup,data.frame(DateReal=as.Date(index(Filtered.Data))))
    days_perFakeYear<-table(year(tmp_DayLookup$DateFake))
    days_perRealYear<-table(year(tmp_DayLookup$DateReal))
    
    if(DataShift=="Yes"){
      Years<-unique(year(tmp_DayLookup$DateFake)) #years in the data
    }else{
      Years<-unique(year(tmp_DayLookup$DateReal)) #years in the data
    }
    
    for(k in 1:length(Years)){ # per year's worth of days for that ROI for that site
      print(paste0(Sites[i],": ROI ",j,"; Year# ",k))
     
      
       if(DataShift=="Yes"){ #get the right datelookup table
        tmp_DayLookup_year<-tmp_DayLookup[year(tmp_DayLookup$DateFake)==Years[k],]
      }else{
        tmp_DayLookup_year<-tmp_DayLookup[year(tmp_DayLookup$DateReal)==Years[k],]  
      }
      
      if(DataShift=="Yes"){ #temporarily trick phenopix by shifting the dates
        Filtered.Data_year<-Filtered.Data[which(year(tmp_DayLookup$DateFake)==Years[k]),]
        index(Filtered.Data_year)<-tmp_DayLookup_year$DateFake[which(year(tmp_DayLookup_year$DateFake)==Years[k])]
      }else{
        Filtered.Data_year<-Filtered.Data[year(index(Filtered.Data))==Years[k],]  
      }
      
      
      
      #to fit all curves and thresholds with no uncertainty estimation,
      exploreFilteredmax <- greenExplore(Filtered.Data_year$max.filtered)
      names(which.min(exploreFilteredmax$rmses))#fitting model with the lowest RMSE 
      plotExplore(exploreFilteredmax)
      assign(paste0("exploreFiltered_",Sites[i],"_ROI",j,"_Year",k),exploreFilteredmax)
      #write to file the explored combinations and uncertainty plots around fit
      jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_Year",k,"_FitPhenoCombos.jpg"))
      plotExplore(exploreFilteredmax)
      dev.off()
      
      
      ## fit with spline with uncertainty
      fitted_withSpline <- SplineFit(Filtered.Data_year$max.filtered, uncert=TRUE, nrep=100)
      fitted_SpKl <- exploreFilteredmax$spline.klosterman# greenProcess(Filtered.Data_year$max.filtered,'spline', 'klosterman',plot=FALSE)
      
      uncertainty.data <- fitted_withSpline$uncertainty$predicted
      sf <- quantile(Filtered.Data_year$max.filtered, na.rm=TRUE, prob=c(0.05, 0.95))
      
      #fit cam's custom trs %15,50, and 90% on either side of peak ; these percentages are based on Bolton et al. 2020
      trs_15<-try(cam_trs_cat(fitted_withSpline$fit$predicted,0.15))
      trs_50<-try(cam_trs_cat(fitted_withSpline$fit$predicted,0.5))
      trs_90<-try(cam_trs_cat(fitted_withSpline$fit$predicted,0.9))
      if(class(trs_15)=="try-error"|class(trs_50)=="try-error"|class(trs_90)=="try-error"){
        trs_pcts<-data.frame(Event=c("SoS_0.15","PoS","EoS_0.15","SoS_0.5","EoS_0.5","SoS_0.9","EoS_0.9"),
                             RealDate=as.Date("2999-01-01"),
                             Value=NA)
        
      }else{
        trs_pcts<-unique(rbind(trs_15,trs_50,trs_90))
    
        if(DataShift=="Yes"){
          trs_pcts$RealDate<-tmp_DayLookup$DateReal[match(trs_pcts$FakeDate,tmp_DayLookup$DateFake)]  
        }else{
          trs_pcts$RealDate<-tmp_DayLookup_year$DateReal[match(trs_pcts$FakeDate,tmp_DayLookup_year$doy) ]  
        }
      
        trs_pcts<-trs_pcts[,c(1,4,3)]
      }
      assign(paste0("FittedSpTRS_",Sites[i],"_ROI",j,"_Year",k),trs_pcts)
      
      thresholds <- NULL
      for (a in 1:dim(uncertainty.data)[2]) { # ie for each replicate of nrep
        ## a passage to correctly format the input data (predicted uncertainty from each rep)
        formatted.data <- list()
        formatted.data$predicted <- uncertainty.data[,a]
        
        ## convert time index into doy
        if(DataShift=="Yes"){
          index(formatted.data$predicted) <- as.numeric(format(index(formatted.data$predicted), '%j'))
        }
        
        ## tmp.column contains thhe four thresholds extracted for a single uncertainty vector 
        tmp.column <- try(suppressWarnings(PhenoKl(x=NULL, fit=formatted.data, uncert=TRUE, breaks=3, sf=sf)))
        if (class(tmp.column)=='try-error') tmp.column <- rep(NA, 4)
        thresholds <- cbind(thresholds, tmp.column)
      }#a
      
      ## threshold is the final matrix, with 4 rows and 100 columns, which you can summarize as follows
      quantiles <- c(0.05, 0.95) ## default as in input to greenProcess
      returned <- as.data.frame(apply(thresholds, 1, function(x) quantile(x, c(quantiles[1], .5, quantiles[2]), na.rm=TRUE)))
      thresholds.t <- t(thresholds)
      names(thresholds.t) <- names(returned)
      rownames(thresholds.t) <- NULL
      fitted_SpKl$metrics<-returned
      fitted_SpKl$uncertainty.df<-as.data.frame(thresholds.t)
      
      
      #reindex results from fitted_SpKl 
      if(DataShift=="Yes"){ #either with the correct dates if data were shifted
        # by the correct dates
        for(l in 1:length(fitted_SpKl$metrics)){ #for greenup, maturity, senecsence and dormancy each,
          fitted_SpKl$uncertainty.df[[l]]<-tmp_DayLookup_year$DateReal[match(fitted_SpKl$uncertainty.df[[l]],tmp_DayLookup_year$doyFake)]
            #tmp_DayLookup$DateReal[match(as.Date(index(Filtered.Data_year)[match(unlist(fitted_SpKl$uncertainty.df[l]),as.numeric(format(index(Filtered.Data_year), "%j")))]),tmp_DayLookup$DateFake)]
          
          for(m in 1:length(fitted_SpKl$metrics[[l]])){ #for each of the quantiles,
            fitted_SpKl$metrics[[l]][m]<-tmp_DayLookup_year$DateReal[match(fitted_SpKl$metrics[[l]][m],tmp_DayLookup_year$doyFake)]
              #tmp_DayLookup$DateReal[which(tmp_DayLookup$DateFake %in%  as.Date(index(Filtered.Data_year)[which(as.numeric(format(index(Filtered.Data_year), "%j"))==unlist(round(fitted_SpKl$metrics[l]))[m])])==1)]
          } #m
          fitted_SpKl$metrics[[l]]<-as.Date( fitted_SpKl$metrics[[l]])
        } #l
        
        index(fitted_SpKl$fit$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fitted_SpKl$fit$fit$predicted),tmp_DayLookup_year$doyFake)]
        index(fitted_SpKl$data)<-index(fitted_SpKl$fit$fit$predicted)
        assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_Year",k),fitted_SpKl)
        
        index(fitted_withSpline$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fitted_withSpline$fit$predicted),tmp_DayLookup_year$DateFake)]
        
      }else{#but if not, then just change of from the DOY to the real date
        
        for(l in 1:length(fitted_SpKl$metrics)){ #for greenup, maturity, senecsence and dormancy each,
          fitted_SpKl$uncertainty.df[[l]]<-tmp_DayLookup_year$DateReal[match(fitted_SpKl$uncertainty.df[[l]],tmp_DayLookup_year$doy)]
          
          for(m in 1:length(fitted_SpKl$metrics[[l]])){ #for each of the quantiles,
            fitted_SpKl$metrics[[l]][m]<-tmp_DayLookup_year$DateReal[match(fitted_SpKl$metrics[[l]][m],tmp_DayLookup_year$doy)]
          }#m
          fitted_SpKl$metrics[[l]]<-as.Date( fitted_SpKl$metrics[[l]])
        } #l
        
        index(fitted_SpKl$fit$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fitted_SpKl$fit$fit$predicted),tmp_DayLookup_year$doy)]
        index(fitted_SpKl$data)<-index(fitted_SpKl$fit$fit$predicted)
        assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_Year",k),fitted_SpKl)
        
        index(fitted_withSpline$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fitted_withSpline$fit$predicted),tmp_DayLookup_year$doy)]
      } #else
      
       
      # fit with klosterman, if possible
      rm(fit.complete) #remove  any old ones
      tryfit<- function(){
        fit.complete<-greenProcess(ts= Filtered.Data_year$max.filtered,
                                   fit = 'klosterman', #'gu',spline', # fitting method 
                                   threshold= 'klosterman',#'gu', #phenophase method
                                   plot = TRUE,
                                   uncert = TRUE,
                                   quantiles=c(0.05,0.95),
                                   nrep = 100)
        
        # reindex the results from fit.complete 
        if(DataShift=="Yes"){ #either by the correct dates
          for(l in 1:length(fit.complete$metrics)){ #for greenup, maturity, senecsence and dormancy each,
            fit.complete$uncertainty.df[[l]]<-tmp_DayLookup_year$DateReal[match(fit.complete$uncertainty.df[[l]],tmp_DayLookup_year$doyFake)]
             
            for(m in 1:length(fit.complete$metrics[[l]])){ #for each of the quantiles,
              fit.complete$metrics[[l]][m]<-tmp_DayLookup_year$DateReal[match(fit.complete$metrics[[l]][m],tmp_DayLookup_year$doyFake)]
            }#m 
            fit.complete$metrics[[l]]<-as.Date( fit.complete$metrics[[l]])
          }#l
         
          index(fit.complete$fit$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fit.complete$fit$fit$predicted),tmp_DayLookup_year$doyFake)]
          index(fit.complete$data)<-index(fit.complete$fit$fit$predicted)
          index(fit.complete$fit$uncertainty$predicted)<-index(fit.complete$fit$fit$predicted)
          
        }else{ # or just with the date rather than doy
          for(l in 1:length(fit.complete$metrics)){ #for greenup, maturity, senecsence and dormancy each,
            fit.complete$uncertainty.df[[l]]<-tmp_DayLookup_year$DateReal[match(fit.complete$uncertainty.df[[l]],tmp_DayLookup_year$doy)]
            
            for(m in 1:length(fit.complete$metrics[[l]])){ #for each of the quantiles,
              fit.complete$metrics[[l]][m]<-tmp_DayLookup_year$DateReal[match(fit.complete$metrics[[l]][m],tmp_DayLookup_year$doy)]
            }#m
            fit.complete$metrics[[l]]<-as.Date(fit.complete$metrics[[l]])
          }#l
          index(fit.complete$fit$fit$predicted)<-tmp_DayLookup_year$DateReal[match(index(fit.complete$fit$fit$predicted),tmp_DayLookup_year$doy)]
          index(fit.complete$data)<-index(fit.complete$fit$fit$predicted)
          index(fit.complete$fit$uncertainty$predicted)<-index(fit.complete$fit$fit$predicted)
          
        }#else
        
             return(fit.complete)
      }
      
      fit.complete<-try(tryfit(),TRUE)
        
      
      #if fit.complete  didnt work,
      if(class(fit.complete)=="try-error"){ 
        fit.complete<-NA
        
        jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_Year",k,"_FittedUncertainty.jpg"))
        par(mfrow=c(3,1))
        plot(index(fitted_withSpline$fit$predicted),fitted_withSpline$fit$predicted,main="Custom trs with spline")
        lines(index(fitted_withSpline$fit$predicted),fitted_withSpline$fit$predicted)
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="SoS_0.15")], col="red", lwd=1, lty=1) #greenup
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="SoS_0.9")], col="green", lwd=1, lty=1) #maturity
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="EoS_0.9")], col="skyblue", lwd=1, lty=1) #senescence
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="EoS_0.15")], col="blue", lwd=1, lty=1) #dormancy
        
        plot(fitted_SpKl)
        plot(1, type="n", xlab="", ylab="",ylim=range(fitted_SpKl$data),xlim=range(index(fitted_SpKl$data)),main="No Klosterman fit due to data-based issue")
        
        dev.off()
      }else{
        #make a plot of the trs, spline and klosterman fits
        assign(paste0("FittedKlKl_",Sites[i],"_ROI",j,"_Year",k),fit.complete)
        jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_Year",k,"_FittedUncertainty.jpg"))
        par(mfrow=c(3,1))
        plot(index(fitted_withSpline$fit$predicted),fitted_withSpline$fit$predicted,main="Custom trs with spline")
        lines(index(fitted_withSpline$fit$predicted),fitted_withSpline$fit$predicted)
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="SoS_0.15")], col="red", lwd=1, lty=1) #greenup
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="SoS_0.9")], col="green", lwd=1, lty=1) #maturity
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="EoS_0.9")], col="blue", lwd=1, lty=1) #senescence
        abline(v = trs_pcts$RealDate[which(trs_pcts$Event=="EoS_0.15")], col="cyan", lwd=1, lty=1) #dormancy
        
        plot(fitted_SpKl)
        plot(fit.complete)
        dev.off()
      }
      
      
    } #k ; per year's worth of days for that ROI for that site
    
    #now combine the "year"-specific phenopix fits (explored and final) for this ROI
    #but only if the data from that "year" is worth it, ie if there's at least half a years worth of days;
    #theyre least likely to be trash model fits
    
    if(DataShift=="Yes"){
      whichYears<-which(days_perFakeYear>(minPropOfYear*365))
    }else{
      whichYears<-which(table(year(tmp_DayLookup$DateReal))>(minPropOfYear*365)) 
    }
     
   
    
    tempNames<-ls(pattern=paste0("exploreFiltered_",Sites[i],"_ROI",j))
    tempKeep<-as.integer(substring(tempNames,nchar(tempNames),nchar(tempNames)))%in% whichYears
    tempNames[tempKeep]
    ROIspecific_ExploredFits_perYear <- do.call("list",mget(tempNames[tempKeep]))
    assign(paste0("exploreFiltered_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_ExploredFits_perYear)
    
    tempNames<-ls(pattern=paste0("FittedKlKl_",Sites[i],"_ROI",j))
    tempKeep<-as.integer(substring(tempNames,nchar(tempNames),nchar(tempNames)))%in%  whichYears
    tempNames[tempKeep]
    if(length(tempNames[tempKeep])>0){
      ROIspecific_FitsKlKl_perYear <- do.call("list",mget(tempNames[tempKeep]))
      assign(paste0("FittedKlKl_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsKlKl_perYear)
    } else{
      ROIspecific_FitsKlKl_perYear <- NA
      assign(paste0("FittedKlKl_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsKlKl_perYear)
    }
    
    tempNames<- ls(pattern=paste0("FittedSpKl_",Sites[i],"_ROI",j))
    tempKeep<-as.integer(substring(tempNames,nchar(tempNames),nchar(tempNames)))%in% whichYears
    tempNames[tempKeep]
    if(length(tempNames[tempKeep])>0){
      ROIspecific_FitsSpKl_perYear <- do.call("list",mget(tempNames[tempKeep]))
      assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsSpKl_perYear)
    } else{
      ROIspecific_FitsSpKl_perYear <- NA
      assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsSpKl_perYear)
    }
    
    tempNames<- ls(pattern=paste0("FittedSpTRS_",Sites[i],"_ROI",j))
    tempKeep<-as.integer(substring(tempNames,nchar(tempNames),nchar(tempNames)))%in% whichYears
    tempNames[tempKeep]
    if(length(tempNames[tempKeep])>0){
      ROIspecific_FitsSpTRS_perYear <- do.call("list",mget(tempNames[tempKeep]))
      assign(paste0("FittedSpTRS_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsSpTRS_perYear)
    } else{
      ROIspecific_FitsSpTRS_perYear <- NA
      assign(paste0("FittedSpTRS_",Sites[i],"_ROI",j,"_perYr"),ROIspecific_FitsSpTRS_perYear)
    }
    
    #reorganize things for those useful/good years in a list called  ROIspecific_Fits_axYears 
    ROIspecific_FitsKlKl_axYears<-list(data=c(),predicted=c(),metrics=c(),uncertaintyPredicted=c(),uncertainty.df=c())
    ROIspecific_FitsSpKl_axYears<-list(data=c(),predicted=c(),metrics=c(),uncertaintyPredicted=c(),uncertainty.df=c())
    ROIspecific_FitsSpTRS_axYears<-list(data=c(),predicted=c(),metrics=c())
    
    if(is.logical(ROIspecific_FitsKlKl_perYear)){ #ie if the KlKl fit didnt work
      ROIspecific_FitsKlKl_axYears <- NA
      assign(paste0("FittedKlKl_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsKlKl_axYears)
      
    }else{
      for(r in 1:length(ROIspecific_FitsKlKl_perYear)){
        print(r)
        ROIspecific_FitsKlKl_axYears$data <- rbind(ROIspecific_FitsKlKl_axYears$data, ROIspecific_FitsKlKl_perYear[[r]]$data)
        ROIspecific_FitsKlKl_axYears$predicted <- rbind(ROIspecific_FitsKlKl_axYears$predicted,ROIspecific_FitsKlKl_perYear[[r]]$fit$fit$predicted)
        ROIspecific_FitsKlKl_axYears$metrics[[r]] <- ROIspecific_FitsKlKl_perYear[[r]]$metrics
        ROIspecific_FitsKlKl_axYears$uncertaintyPredicted[[r]]<-ROIspecific_FitsKlKl_perYear[[r]]$fit$uncertainty$predicted
        ROIspecific_FitsKlKl_axYears$uncertainty.df[[r]] <- ROIspecific_FitsKlKl_perYear[[r]]$uncertainty.df
      }
      assign(paste0("FittedKlKl_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsKlKl_axYears)
    }
    
    if(is.logical(ROIspecific_FitsSpKl_perYear)){ #ie if this didnt work
      ROIspecific_FitsSpKl_axYears <- NA
      assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsSpKl_axYears)
      
    }else{
      for(r in 1:length(ROIspecific_FitsSpKl_perYear)){ 
        ROIspecific_FitsSpKl_axYears$data <- rbind(ROIspecific_FitsSpKl_axYears$data, ROIspecific_FitsSpKl_perYear[[r]]$data )
        ROIspecific_FitsSpKl_axYears$predicted <- rbind(ROIspecific_FitsSpKl_axYears$predicted,ROIspecific_FitsSpKl_perYear[[r]]$fit$fit$predicted)
        ROIspecific_FitsSpKl_axYears$metrics[[r]] <- ROIspecific_FitsSpKl_perYear[[r]]$metrics
        ROIspecific_FitsSpKl_axYears$uncertaintyPredicted[[r]]<-ROIspecific_FitsSpKl_perYear[[r]]$fit$uncertainty$predicted
        ROIspecific_FitsSpKl_axYears$uncertainty.df[[r]] <- ROIspecific_FitsSpKl_perYear[[r]]$uncertainty.df
      }
      assign(paste0("FittedSpKl_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsSpKl_axYears)
    }
    
    if(is.logical(ROIspecific_FitsSpTRS_perYear)){ #ie if this didnt work
      ROIspecific_FitsSpTRS_axYears <- NA
      assign(paste0("FittedSpTRS_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsSpTRS_axYears)
      
    }else{
      for(r in 1:length(ROIspecific_FitsSpTRS_perYear)){ 
        ROIspecific_FitsSpTRS_axYears$data <- rbind(ROIspecific_FitsSpTRS_axYears$data, ROIspecific_FitsSpKl_perYear[[r]]$data )
        ROIspecific_FitsSpTRS_axYears$predicted <- rbind(ROIspecific_FitsSpTRS_axYears$predicted,ROIspecific_FitsSpKl_perYear[[r]]$fit$fit$predicted)
        ROIspecific_FitsSpTRS_axYears$metrics[[r]] <- ROIspecific_FitsSpTRS_perYear[[r]]
      }
      assign(paste0("FittedSpTRS_",Sites[i],"_ROI",j,"_axYrs"),ROIspecific_FitsSpTRS_axYears)
    }
    
    ##prep for plots
    whichWorked<-1+c(is.logical(ROIspecific_FitsKlKl_axYears),is.logical(ROIspecific_FitsSpKl_axYears),
                     is.logical(ROIspecific_FitsSpTRS_axYears))
    # if  none of the fitting attempts worked
    if(sum(whichWorked)==0){
      gPlot1<-"No SpTRS was successfull for this ROI in this Year at this Site"
      gPlot2<-"No SpKl was successfull for this ROI in this Year at this Site"
      gPlot3<-"No KlKl was successfull for this ROI in this Year at this Site"
      assign(paste0("Plot_Fitted_KlKl","_ROI",j),gPlot3)
      assign(paste0("Plot_Fitted_SpKl","_ROI",j),gPlot2)
      assign(paste0("Plot_Fitted_SpTRS","_ROI",j),gPlot1)
    }
    tmp_plot_pred<-list()
    tmp_plot_dat<-list()
    tmp_plot_uncert<-list()
    tmp_plot_metrics<-list()
    
    #if SpKl worked
    if(whichWorked[2]==1){ 
      tmp_plot_pred$SpKl<-fortify.zoo(ROIspecific_FitsSpKl_axYears$predicted,melt=TRUE)
      tmp_plot_pred$SpKl<-cbind(rep(Sites[i],nrow(tmp_plot_pred$SpKl)),rep(j,nrow(tmp_plot_pred$SpKl)),tmp_plot_pred$SpKl[-2])
      colnames(tmp_plot_pred$SpKl)<-c("Site","ROI","Date","Pred")
      
      #prep the data to plot
      tmp_plot_dat$SpKl<-fortify.zoo(ROIspecific_FitsSpKl_axYears$data,melt=TRUE)
      tmp_plot_dat$SpKl<-cbind(rep(Sites[i],nrow(tmp_plot_dat$SpKl)),rep(j,nrow(tmp_plot_dat$SpKl)),tmp_plot_dat$SpKl[-2])
      colnames(tmp_plot_dat$SpKl)<-c("Site","ROI","Date","Value")
      
      #prep the uncertainty to plot
      tmp_plot_uncert$SpKl<-fortify.zoo(ROIspecific_FitsSpKl_axYears$uncertaintyPredicted,melt=TRUE)
      tmp_plot_uncert$SpKl<-cbind(rep(Sites[i],nrow(tmp_plot_uncert$SpKl)),rep(j,nrow(tmp_plot_uncert$SpKl)),tmp_plot_uncert$SpKl[-2])
      
      #prep the metrics to plot
      tmp_plot_metrics$SpKl<-transpose(do.call("cbind",ROIspecific_FitsSpKl_axYears$metrics))
      
      tmp_plot_metrics$SpKl <-as.data.frame(lapply(tmp_plot_metrics$SpKl,as.Date,origin="1970-01-01"))
      tmp_plot_metrics$SpKl<- cbind(rep(Sites[i],length(ROIspecific_FitsSpKl_perYear)),rep(j,length(ROIspecific_FitsSpKl_perYear)),
                                    rep(c("Greenup","Maturity","Senescence","Dormancy"),length(ROIspecific_FitsSpKl_perYear)),tmp_plot_metrics$SpKl)
      rownames(tmp_plot_metrics$SpKl)<-unlist(lapply(seq(1,length(ROIspecific_FitsSpKl_perYear),by=1),paste0,c("Greenup","Maturity","Senescence","Dormancy")))
      colnames(tmp_plot_metrics$SpKl)<-c("Site","ROI","Event","5%","50%","95%")
    }else{
      jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_allYrs_Plot_Fitted_SpKl.jpg"))
      plot(1, type="n", xlab="", ylab="",ylim=range(0,1),xlim=range(0,1),main="No Spline fit due to data-based issues")
      dev.off()
    }
    
    #if KlKl worked
    if(whichWorked[1]==1){
      #prep the fit to plot
      tmp_plot_pred$KlKl<-fortify.zoo(ROIspecific_FitsKlKl_axYears$predicted,melt=TRUE)
      tmp_plot_pred$KlKl<-cbind(rep(Sites[i],nrow(tmp_plot_pred$KlKl)),rep(j,nrow(tmp_plot_pred$KlKl)),tmp_plot_pred$KlKl[-2])
      colnames(tmp_plot_pred$KlKl)<-c("Site","ROI","Date","Pred")
      #assign(paste0("Fitted_preds_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_pred)
      
      #prep the data to plot
      tmp_plot_dat$KlKl<-fortify.zoo(ROIspecific_FitsKlKl_axYears$data,melt=TRUE)
      tmp_plot_dat$KlKl<-cbind(rep(Sites[i],nrow(tmp_plot_dat$KlKl)),rep(j,nrow(tmp_plot_dat$KlKl)),tmp_plot_dat$KlKl[-2])
      colnames(tmp_plot_dat$KlKl)<-c("Site","ROI","Date","Value")
      #assign(paste0("Plotted_dat_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_dat)
      
      #prep the uncertainty to plot
      tmp_plot_uncert$KlKl<-fortify.zoo(do.call("rbind",ROIspecific_FitsKlKl_axYears$uncertaintyPredicted),melt=TRUE)
      tmp_plot_uncert$KlKl<-cbind(rep(Sites[i],nrow(tmp_plot_uncert$KlKl)),rep(j,nrow(tmp_plot_uncert$KlKl)),tmp_plot_uncert$KlKl[-2])
      colnames(tmp_plot_uncert$KlKl)<-c("Site","ROI","Date","Value")
      #assign(paste0("Fitted_uncertainty_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_uncert)
      
      #prep the metrics to plot
      tmp_plot_metrics$KlKl<-transpose(do.call("cbind",ROIspecific_FitsKlKl_axYears$metrics))
      tmp_plot_metrics$KlKl <-as.data.frame(lapply(tmp_plot_metrics$KlKl,as.Date,origin="1970-01-01"))
      tmp_plot_metrics$KlKl<- cbind(rep(Sites[i],length(ROIspecific_FitsKlKl_perYear)),rep(j,length(ROIspecific_FitsKlKl_perYear)),
                                    rep(c("Greenup","Maturity","Senescence","Dormancy"),length(ROIspecific_FitsKlKl_perYear)),tmp_plot_metrics$KlKl)
      rownames(tmp_plot_metrics$KlKl)<-unlist(lapply(seq(1,length(ROIspecific_FitsKlKl_perYear),by=1),paste0,c("Greenup","Maturity","Senescence","Dormancy")))
      
      colnames(tmp_plot_metrics$KlKl)<-c("Site","ROI","Event","5%","50%","95%")
      #assign(paste0("Fitted_metrics_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_metrics)
      
    }else{
      jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_allYrs_Plot_Fitted_KlKl.jpg"))
      plot(1, type="n", xlab="", ylab="",ylim=range(0,1),xlim=range(0,1),main="No Klosterman fit due to data-based issues")
      dev.off()
    }
    #if SpTRS worked
    if(whichWorked[3]==1){
      tmp_plot_pred$SpTRS<-fortify.zoo(ROIspecific_FitsSpKl_axYears$predicted,melt=TRUE)
      tmp_plot_pred$SpTRS<-cbind(rep(Sites[i],nrow(tmp_plot_pred$SpTRS)),rep(j,nrow(tmp_plot_pred$SpTRS)),tmp_plot_pred$SpTRS[-2])
      colnames(tmp_plot_pred$SpTRS)<-c("Site","ROI","Date","Pred")
      
      #prep the data to plot
      tmp_plot_dat$SpTRS<-fortify.zoo(ROIspecific_FitsSpKl_axYears$data,melt=TRUE)
      tmp_plot_dat$SpTRS<-cbind(rep(Sites[i],nrow(tmp_plot_dat$SpTRS)),rep(j,nrow(tmp_plot_dat$SpTRS)),tmp_plot_dat$SpTRS[-2])
      colnames(tmp_plot_dat$SpTRS)<-c("Site","ROI","Date","Value")
      
      #prep the uncertainty to plot
      tmp_plot_uncert$SpTRS<-NA
      
      #prep the metrics to plot
      tmp_plot_metrics$SpTRS<-do.call("rbind",ROIspecific_FitsSpTRS_axYears$metrics)[,1:2]
      tmp_plot_metrics$SpTRS<-tmp_plot_metrics$SpTRS[which(tmp_plot_metrics$SpTRS$Event!="PoS"&
                                                             tmp_plot_metrics$SpTRS$Event!="SoS_0.5"&
                                                             tmp_plot_metrics$SpTRS$Event!="EoS_0.5"),]
      tmp_plot_metrics$SpTRS$Event<-as.character(tmp_plot_metrics$SpTRS$Event)
      tmp_plot_metrics$SpTRS$Event[which(tmp_plot_metrics$SpTRS$Event=="SoS_0.15")]<-"Greenup"
      tmp_plot_metrics$SpTRS$Event[which(tmp_plot_metrics$SpTRS$Event=="EoS_0.15")]<-"Dormancy"
      tmp_plot_metrics$SpTRS$Event[which(tmp_plot_metrics$SpTRS$Event=="SoS_0.9")]<-"Maturity"
      tmp_plot_metrics$SpTRS$Event[which(tmp_plot_metrics$SpTRS$Event=="EoS_0.9")]<-"Senescence"
      tmp_plot_metrics$SpTRS$Event<-as.factor(tmp_plot_metrics$SpTRS$Event)
      tmp_plot_metrics$SpTRS<- cbind(rep(Sites[i],nrow(tmp_plot_metrics$SpTRS)),rep(j,nrow(tmp_plot_metrics$SpTRS)),
                                     tmp_plot_metrics$SpTRS)
      colnames(tmp_plot_metrics$SpTRS)<-c("Site","ROI","Event","50%")
      tmp_plot_metrics$SpTRS<- tmp_plot_metrics$SpTRS[order(tmp_plot_metrics$SpTRS$`50%`),]
      #rownames(tmp_plot_metrics$SpTRS)<-unlist(lapply(seq(1,length(ROIspecific_FitsSpKl_perYear),by=1),paste0,c("Greenup","Maturity","Senescence","Dormancy")))
      
    }else{
      jpeg(paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_allYrs_Plot_Fitted_SpTRS.jpg"))
      plot(1, type="n", xlab="", ylab="",ylim=range(0,1),xlim=range(0,1),main="No Spline fit due to data-based issues")
      dev.off()
    }
    
    assign(paste0("Fitted_preds_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_pred)
    assign(paste0("Plotted_dat_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_dat)
    assign(paste0("Fitted_uncertainty_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_uncert)
    assign(paste0("Fitted_metrics_",Sites[i],"_ROI",j,"_allYrs"),tmp_plot_metrics)
    
    
    
    
    #distinguish between datapoints used and those discarded; could prob also just add a column in Filtered.Data saying as much.
    if(DataShift=="Yes"){
      data_modeled<-Filtered.Data[which(year(tmp_DayLookup$DateFake) %in% as.numeric(attributes(which(days_perFakeYear>(minPropOfYear*365)))$names)),]
      `%notin%` <- Negate(`%in%`)
      data_notModeled<-Filtered.Data[which(year(tmp_DayLookup$DateFake) %notin% as.numeric(attributes(which(days_perFakeYear>(minPropOfYear*365)))$names)),]
       
    }else{
      
      data_modeled<-Filtered.Data[which(year(tmp_DayLookup$DateReal) %in% as.numeric(attributes(which(days_perRealYear>(minPropOfYear*365)))$names)),]
      `%notin%` <- Negate(`%in%`)
      data_notModeled<-Filtered.Data[which(year(tmp_DayLookup$DateReal) %notin% as.numeric(attributes(which(days_perRealYear>(minPropOfYear*365)))$names)),]
      
    }
    ROInames<-c(1:nroi,1)
    
    
    #create and save plots
    for(method in 1:length(tmp_plot_metrics)){ # for KlKl,  SpKl, SpTRS separately
      metric_colors<-rep(c("red4","darkgreen","dodgerblue4","steelblue4"),nrow(tmp_plot_metrics[[method]])/4)
      metric_bgcolors<-rep(c("red","green","dodgerblue1","steelblue1"),nrow(tmp_plot_metrics[[method]])/4)
      
      gPlot<-ggplot()+
        xlim(low=min(as.Date(index(Filtered.Data))), high=max(as.Date(index(Filtered.Data))))+
        #ylim(low=0.3, high=0.6)+
        ylim(low=min(tmp_plot_dat[[method]]$Value)-0.01, high=max(tmp_plot_dat[[method]]$Value)+0.1)
      if(names(tmp_plot_metrics[method])!="SpTRS"){
        gPlot<-gPlot+geom_rect(data=tmp_plot_metrics[[method]],aes(xmin = `5%`, xmax = `95%`, ymin = -Inf, ymax = Inf),alpha=0.9,fill=metric_bgcolors)#,fill = "red",color = "red4")+
      }
      if(names(tmp_plot_metrics[method])=="KlKl"){
        gPlot<-gPlot+geom_line(data=tmp_plot_uncert[[method]],aes(x=Date,y=Value),col="gray")
      }
      gPlot<-gPlot+
        geom_line(data=tmp_plot_pred[[method]],aes(x=Date,y=Pred))+
        geom_vline(data=tmp_plot_metrics[[method]],aes(xintercept = `50%`),color=metric_colors)+
        geom_point(data=data_modeled,aes(x=as.Date(index(data_modeled)),y=max.filtered))+
        scale_shape_identity()+
        geom_text(data = tmp_plot_metrics[[method]], mapping = aes(x=`50%`,label = `50%`, y = max(tmp_plot_dat[[method]]$Value)+0.01), angle = 90, hjust = 0)+
        ggtitle(paste("ROI",ROInames[j],sep=" ")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.title.x = element_blank(),axis.title.y = element_blank())
      if(dim(data_notModeled)[1]>0){
        gPlot<-gPlot+geom_point(data=data_notModeled,aes(x=as.Date(index(data_notModeled)),y=max.filtered,shape=1))
      }
      
      #gPlot
      ggsave(gPlot, width=15,height=10,units="in" ,
             file=paste0("./FitAndPheno_plots/",Sites[i],"/",Sites[i],"_ROI",j,"_allYrs_Plot_Fitted_",names(tmp_plot_metrics)[method],".jpg"))
      assign(paste0("Plot_Fitted_",names(tmp_plot_metrics)[method],"_ROI",j),gPlot)
      
    }
    
    
    
  }# j ; per roi for that site
  #Into a single list, combine metrics (mean and quantiles of phenodates) from the Spline_Klosterman and Klosterman_Klosterman 
  
  AllMetrics_Fitted_AllROIs_overTime <- sapply(c("KlKl","SpKl","SpTRS"),function(x) NULL)
  for(s in 1:length(ls(pattern="Fitted_metrics_"))){ #for each ROI
    temp<-get(ls(pattern="Fitted_metrics_")[s])
    names(temp)
    temp
    if("SpKl"%in% names(temp)){
      AllMetrics_Fitted_AllROIs_overTime$SpKl<-rbind(AllMetrics_Fitted_AllROIs_overTime$SpKl, temp$SpKl)
      rownames(AllMetrics_Fitted_AllROIs_overTime$SpKl)<-NULL
    }
    
    if("KlKl"%in% names(temp)){
      AllMetrics_Fitted_AllROIs_overTime$KlKl<-rbind( AllMetrics_Fitted_AllROIs_overTime$KlKl,temp$KlKl)
      rownames(AllMetrics_Fitted_AllROIs_overTime$KlKl)<-NULL
    }
    
    if("SpTRS"%in% names(temp)){
      AllMetrics_Fitted_AllROIs_overTime$SpTRS<-rbind(AllMetrics_Fitted_AllROIs_overTime$SpTRS,temp$SpTRS)
      rownames(AllMetrics_Fitted_AllROIs_overTime$SpTRS)<-NULL
    }
  }
  
  
  plist_KlKl<- do.call("list",mget(grep("Plot_Fitted_KlKl",names(.GlobalEnv),value=TRUE)))
  plist_SpKl<- do.call("list",mget(grep("Plot_Fitted_SpKl",names(.GlobalEnv),value=TRUE)))
  plist_SpTRS<- do.call("list",mget(grep("Plot_Fitted_SpTRS",names(.GlobalEnv),value=TRUE)))
  plist_all<- do.call("list",mget(grep("Plot_Fitted_",names(.GlobalEnv),value=TRUE)))
  
  gPlot_AllROIs_KlKl_overTime<- grid.arrange(grobs=plist_KlKl,
                                             ncol=1,
                                             top =textGrob(Sites[i], gp=gpar(fontsize=20)) ,
                                             bottom=textGrob("Date", gp=gpar(fontsize=18)),
                                             left=textGrob("Filtered Relative Green Chromatic Coordinates",gp=gpar(fontsize=18),rot=90))
  gPlot_AllROIs_SpKl_overTime<- grid.arrange(grobs=plist_SpKl,
                                             ncol=1,
                                             top =textGrob(Sites[i], gp=gpar(fontsize=20)) ,
                                             bottom=textGrob("Date", gp=gpar(fontsize=18)),
                                             left=textGrob("Filtered Relative Green Chromatic Coordinates",gp=gpar(fontsize=18),rot=90))
  gPlot_AllROIs_SpTRS_overTime<- grid.arrange(grobs=plist_SpTRS,
                                              ncol=1,
                                              top =textGrob(Sites[i], gp=gpar(fontsize=20)) ,
                                              bottom=textGrob("Date", gp=gpar(fontsize=18)),
                                              left=textGrob("Filtered Relative Green Chromatic Coordinates",gp=gpar(fontsize=18),rot=90))
  
  
  ggsave(gPlot_AllROIs_KlKl_overTime,width=8,height=12,units="in" ,file=paste0("./FitAndPheno_plots/",Sites[i],"/","Plot_Fitted_KlKl_",Sites[i],"_allYrs.jpg"))
  ggsave(gPlot_AllROIs_SpKl_overTime,width=8,height=12, file=paste0("./FitAndPheno_plots/",Sites[i],"/","Plot_Fitted_SpKl_",Sites[i],"_allYrs.jpg"))
  ggsave(gPlot_AllROIs_SpTRS_overTime,width=8,height=12, file=paste0("./FitAndPheno_plots/",Sites[i],"/","Plot_Fitted_SpTRS_",Sites[i],"_allYrs.jpg"))
  
  
  #save the explored and final Fits across years for all rois
  #save(thingsToSaveToFile,file = paste0(OutputDir,"Rdata/","PhenopixAnalysisResults_",Sites[i],".RData"))
  
  save.image(file = paste0("./Rdata/","PhenopixAnalysisResults_",Sites[i],".RData"))
  # rm(list=setdiff(ls(), c("OutputDir",
  #                         "Sites","DataShift",
  #                         "i","cam_trs_cat","minPropOfYear")))
}# per site 


