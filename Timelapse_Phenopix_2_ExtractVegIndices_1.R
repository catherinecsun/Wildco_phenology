##############################################################
# Timelapse_Phenopix_2_ExtractVegIndices.R
# Extracts the Veg Indices per site per season
# option available to do the phenophase analysis per season
# but may be extra work if ultimately you want to analyse all deployments together
# an if() function is used to control whether those phenophase analyses per deployment are conducted
##############################################################

####--- load packages ####
library(phenopix)
library(jpeg)
library(rasterImage)
library(stringr)

###---- Define and Set things THIS NEEDS TO BE SPECIFIC TO YOUR PROJECT AND YOUR FOLDER STRUCTURE. ####
OutputDir <- "./Output/"
locOfImgs <- "./Timelapse_photos/" #the big folder where all the timelapse photos are (with all sites, all seasons)

#define the site(s) you want to do the analysis on; you can put in mutiple, but they can take a while so 
# probably best to just do one first, and then ramp up to a handful at a time to  better keep an eye on things
Site <-c("CATH08")


#define how many of the ROIs you want to do
# i have also moved this to inside the loop below because the nroi could change for each site
#each site may hae a different number of ROIs drawn. so the code automatically figures that out now
nroi_intent<-1


###---Loop through to do Veg Extraction ####
# for each site
# pull out the Veg Indices and other info: ie, the raw + RGB digital numbers, chromatic coordinates
# Should be noted that for date.code argument in extractVIs()
# that the naming convention of files must be set up properly for the function to work. Refer to vignette.
for(a in 1:length(Site)){
  print(Site[a])
  setwd(paste0(OutputDir,"ROI/",Site[a],"/"))
  
  
  #figure out how many ROIs the analyses need to be done for
  site_nameLength<-nchar(strsplit(list.files(pattern = ".RData"),"_")[[1]][1])
  nroi<-as.numeric(substring(list.files(pattern = ".RData"),site_nameLength+2,site_nameLength+2))
  
  #print(paste0("# ROIs: ",nroi))
  # if(nroi_intent<nroi){
  #   print(paste0("there seem to be ", nroi, "ROI(s),but we are only doing" , nroi_intent))
  #   nroi_intent<-nroi
  # }
 
  #figure out how many deployment periods the analyses need to be done for
  # i.e., for each site define the periods you want to do the veg analysis on 
 
  tmp_whichToExtract<-lapply(strsplit(list.files(pattern = ".RData"),"_"),length)[[1]]
  SeasonPerSite<-sapply( strsplit(list.files(pattern = ".RData"),"_"), "[", tmp_whichToExtract )
  SeasonPerSite<- substring(SeasonPerSite,1,nchar(SeasonPerSite)-6) #get rid of the ".RData" at the end
  
  #inelegantly go back up 3 levels to the big Folder
  setwd('..')
  setwd('..')
  setwd('..')

  #create a site-specific folder in the VI folder to put the site specific VI results in
  dir.create(paste0(OutputDir,"VI/",Site[a]))
  
  #now, loop through the analysis for each season/deployment period
for(i in 1:length(SeasonPerSite)){
  print(SeasonPerSite[i])
  Season<-SeasonPerSite[i]
  
  #create a temp copy of the ROI.rdata with a generic name
  trycopy<-try(file.copy(paste0(OutputDir,"ROI/",Site[a] ,"/",Site[a],"_",nroi,"ROI","_",Season,".RData"), 
            paste0(OutputDir,"ROI/",Site[a],"/roi.data.RData"),overwrite = TRUE))  
  while (trycopy==0) { #if that didn't work, go on to the next season/deployment
    print("didnt work")
    i <- i+1
    Season<-SeasonPerSite[i]
    print(SeasonPerSite[i])
    trycopy<-try(file.copy(paste0(OutputDir,"ROI/",Site[a] ,"/",Site[a],"_",nroi,"ROI","_",Season,".RData"), 
                           paste0(OutputDir,"ROI/",Site[a],"/roi.data.RData"),overwrite = TRUE)) #create a temp copy of the ROI.rdata with the generic name 
    
  }
  
  load(paste(paste0(OutputDir,"ROI/",Site[a],"/"), "/roi.data.Rdata", sep = ""))
  length(names(roi.data))
  if(nroi<length(names(roi.data))){
    print(paste0("there are ", length(names(roi.data)), " ROI(s),but we are only doing the " , nroi,"nd one"))
    nroi_intent<-nroi
  }
  ###---extract the VegIndices stuff####
  #this extractVIs can take a while
  VI.data <-extractVIs(paste0(locOfImgs,Site[a],"/",Season,"/"), # image path
                       paste0(OutputDir,"ROI/",Site[a],"/"), # roi path to get the stored roi 
                       vi.path = paste0(OutputDir,"VI/",Site[a],"/"), # where VI will be saved
                       roi.name =  paste0(Site[a],"_",Season,"_roi",seq(1,nroi,1)), #name of the roi in the ROI.rdata
                       plot = TRUE, begin = NULL, spatial = FALSE, date.code= 'yyyy-mm-dd__HH-MM', npixels = 1,
                       file.type = ".jpg", bind = FALSE)
  
  
  #There should be one list for each ROI, and check that  proper columns (18) are included
  #length(VI.data)
  #names(VI.data)
  #View(VI.data[[1]])
  #names(VI.data[[1]])
  
  # Turn VI.data from a list into data frame and save to file
  VI.data.df <- do.call("rbind",VI.data)
  seq(1,length(VI.data),1)
  ROI<-c()
  for(z in 1:length(VI.data)){
    temp<-rep(z,nrow(VI.data[[z]]))
    ROI<-c(ROI,temp)
  }
  VI.data.df <-cbind(ROI,VI.data.df)
  assign(paste0("VI.data.df_",Site[a],"_",Season),VI.data.df)
  
  #rename stuff 
   #first VI plot
  file.rename(from=paste0(OutputDir,"VI/", Site[a],"/",Site[a],"_",Season,"_roi1_roi_VI_plot",".png"),
              to=paste0(OutputDir,"VI/", Site[a],"/", Site[a],"_roi1","_VIplot_",Season,".png"))
  #second VI plot: if there is more than 1 roi, I cant figure out how to give the VI plots unique names. so rename one ones after first from NA to 
  if(nroi>1){
    for(q in 2:nroi){
    file.rename(from=paste0(OutputDir,"VI/", Site[a],"/",Site[a],"_",Season,"_roi",q,"_roi_VI_plot",".png"),
                to=paste0(OutputDir,"VI/", Site[a],"/",Site[a],"_roi",q,"_VIplot_",Season,".png"))
    }
  }
  
  #rename the Rdata from generic VI.data.RData to a site-specific name
  file.rename(from=paste0(OutputDir,"VI/",Site[a],"/VI.data.RData"),
              to=paste0(OutputDir,"VI/",Site[a],"/",Site[a],"_",Season,"_VIdata.RData"))
  #remove the generic ROI rdata file we had to temporarily create
  file.remove(paste0(OutputDir,"ROI/",Site[a],"/roi.data.RData"))

 } #analysis loop per deployment
  
  #save the VI data  across all  ROIs and all seasons/deployment periods in a single csv.
  #ls(pattern=paste0("VI.data.df_",Site[a],"_"))
  VIs_perSeason<- do.call("list",mget(ls(pattern=paste0("VI.data.df_",Site[a],"_"))))
  VIs_axSeason<-do.call("rbind",VIs_perSeason)
  VIs_axSeason<-VIs_axSeason[order(VIs_axSeason$date),]
  VIs_axSeason<-VIs_axSeason[order(VIs_axSeason$ROI),]
  assign(paste0(Site[a],"_VI1and2_axSeasons"),VIs_axSeason)  
  write.csv(VIs_axSeason,paste0(OutputDir,"VI/",Site[a],"/", Site[a],"_VI_axSeasons",".csv"))
  
  #create a plot with all the VI's - all ROIs, all seasons
  #NOTE THAT IF THERE ARE MULTIPLE ROIS THEN PLEASE UN-HASHTAG THE CORRESPONDING LINES BELOW SO THAT DATA
  # FROM THOSE ROIS ALSO SHOW UP IN THE PLOT. 
  png(filename = paste0(OutputDir,"VI/",Site[a],"/",Site[a],"_VI_axSeasons_ROIplot",".png"), width = 800, height = 5 * 400, pointsize = 30)
  par(mfrow = c(5, 1))
  par(mar = c(3, 4, 2, 0.5))
  plot(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$r.av[which(VIs_axSeason$ROI==1)], col = adjustcolor("darkred",alpha=1),
       pch = 20, xlab = "", ylab = "R-G-B",xaxt="n", ylim=c(min(c(VIs_axSeason$r.av,VIs_axSeason$g.av,VIs_axSeason$b.av)),
                                                            max(c(VIs_axSeason$r.av,VIs_axSeason$g.av,VIs_axSeason$b.av))),
       main = paste(Site[a],": Color indices in ROIs over time ", sep = ""),
       bty='L' )
  axis.POSIXct(1, at = seq(VIs_axSeason$date[1], VIs_axSeason$date[length(VIs_axSeason$date)], by = "months"), format = "%b")
  points(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$g.av[which(VIs_axSeason$ROI==1)], col = adjustcolor("darkgreen",alpha=1), pch = 20)
  points(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$b.av[which(VIs_axSeason$ROI==1)], col = adjustcolor("darkblue",alpha=1),pch = 20)
 
  # to include more rois in this plot, could come up with a for loop..
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$r.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("red", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$g.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("darkgreen", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$b.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("darkblue", alpha=0.5))
  # 
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$r.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("orangered", alpha=0.7), pch = 20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$g.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("palegreen3", alpha=0.4), pch = 20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$b.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("skyblue", alpha=0.5), pch = 20)
  # 
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$r.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("deeppink", alpha=0.7), pch = 20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$g.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("olivedrab1", alpha=0.7), pch = 20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$b.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("turquoise1", alpha=0.5), pch = 20)
  
   par(xpd=TRUE)
  legend("topleft", legend=c("Entire ROI ","Centered ROI","Left ROI","Right ROI") , pch=c(20,1,20,20), 
         col=c("black","black","gainsboro","gray"), horiz=TRUE,bty="n")
  
  par(mar = c(3, 4, 0.5, 0.5))
  plot(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$ri.av[which(VIs_axSeason$ROI==1)], col = "darkred", 
       ylim=c(min(VIs_axSeason$ri.av),max(VIs_axSeason$ri.av)), pch = 20, xlab = "", ylab = "RI")
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$ri.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("red", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$ri.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("orangered", alpha=0.7),pch=20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$ri.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("deeppink", alpha=0.7),pch=20)
  
  par(mar = c(3, 4, 0.5, 0.5))
  plot(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$gi.av[which(VIs_axSeason$ROI==1)], col = "darkgreen",
       ylim=c(min(VIs_axSeason$gi.av),max(VIs_axSeason$gi.av)), pch = 20, xlab = "", ylab = "GI")
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$gi.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("darkgreen", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$gi.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("palegreen3", alpha=0.4),pch=20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$gi.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("olivedrab1", alpha=0.7),pch=20)
  
  par(mar = c(3, 4, 0.5, 0.5))
  plot(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$bi.av[which(VIs_axSeason$ROI==1)], col = "darkblue", 
       ylim=c(min(VIs_axSeason$bi.av),max(VIs_axSeason$bi.av)), pch = 20, xlab = "", ylab = "BI")
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$bi.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("darkblue", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$bi.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("skyblue", alpha=0.5),pch=20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$bi.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("turquoise1", alpha=0.5),pch=20)
  
  par(mar = c(4, 4, 0.5, 0.5))
  plot(VIs_axSeason$date[which(VIs_axSeason$ROI==1)], VIs_axSeason$bri.av[which(VIs_axSeason$ROI==1)], col = "black",
       ylim=c(min(VIs_axSeason$bri.av),max(VIs_axSeason$bri.av)), pch = 20, xlab = "doy", ylab = "BRI")
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==2)], VIs_axSeason$bri.av[which(VIs_axSeason$ROI==2)], col = adjustcolor("black", alpha=0.5))
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==3)], VIs_axSeason$bri.av[which(VIs_axSeason$ROI==3)], col = adjustcolor("gray72", alpha=1),pch=20)
  # points(VIs_axSeason$date[which(VIs_axSeason$ROI==4)], VIs_axSeason$bri.av[which(VIs_axSeason$ROI==4)], col = adjustcolor("gray47", alpha=1),pch=20)
  
  dev.off()
}#analysis loop per site

