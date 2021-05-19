##############################################################
# Timelapse_Phenopix_1_CreateROIs_Cat.R
# Creates Folders and the ROIs per site per Dplymnt/deployment period
# all photos for a site for the Dplymnt/deploment  you want to analyse should be in the same folder
# name for that folder should be 'Site_Dplymnt' e.g., Algar01_AprToNov2018
# have a big output folder where phenopix output will go.
# it would be really hard to do a single ROI for the entire year/ie multiple Dplymnts becaues the camera moves and the area 
# of view changes so the same exact ROI cant be applied... see <http://graphics.cs.cmu.edu/projects/crossDomainMatching/abhinav-sa11.pdf>
# to get patterns across multiple Dplymnts are desired, we will combine them afterwards see the next script in this sequence.
##############################################################

####--- load packages ####
library(phenopix)
library(jpeg)
library(rasterImage)
library(stringr)

###---- Define and Set things  THIS NEEDS TO BE SPECIFIC TO YOUR PROJECT AND YOUR FOLDER STRUCTURE. ####
# point to where the timelapse images are
# create an Output directory for everything else to go, and point to it
locOfImgs <- "./Timelapse_photos/" 
dir.create("./Output/")
OutputDir <- "./Output/"

### THIS NEEDS TO BE CHANGED BY THE USER
# define Site for which you want to do the analysis/draw the ROI
#do ths one site at a time. 
Site<-"CATH08"

# define the Deploymnts you want to focus on, by looking at the folders
#this will take a little bit of time
#list.dirs(locOfImgs, recursive=FALSE) 
#setwd(OutputDir)
dirs <- grep(Site, list.dirs(recursive=TRUE), value = TRUE)
#Dplymnt <-substring(dirs,105,nchar(dirs)[1])

#get the deployments
Dplymnt <-  sapply( strsplit(dirs,"/"), "[", 4 )
Dplymnt<-Dplymnt[is.na(Dplymnt)==FALSE]

### THIS NEEDS TO BE INPUT BY THE USER
# give a vector the name of the imgs you want to use as references;one img per Dplymnt
#and give them in the order of the Dplymnts
useForROI<-c('CATH08__2019-07-05__12-00-00.JPG')

length(useForROI)==length(Dplymnt)# make sure you have an image per folder/Dplymnt"

### THIS NEEDS TO BE INPUT BY THE USER
# define number of ROIs per image you'll be creating. ONE OR TWO. 
nroi<-1

#will you want to look at the consistency of your ROI placement over time?
# if yes, the code will assume you are drawing quadrilaterals and wont work otherwise. 
# if you want to draw weirdly shaped polygons, you will need to make the code accommodate that shape 
# or add flexibility to allow whatever shape you draw per image... just probably just stick to a quadrilateral.
ROIcheck<-"N" # "N"

#create an empty list in which to put the 4 ROI vertices; so if you want some other polygon, you will have to change/add number of vertices
if(ROIcheck=="Y"){
  ROIvertices<-data.frame(Dplymnt=character(),ROI=integer(),x1=double(),x2=double(),x3=double(),x4=double(),y1=double(),y2=double(),y3=double(),y4=double())
  ROIvertices$Dplymnt<-as.character(ROIvertices$Dplymnt)
}

# determine if necessary output folders exist, and create them if not
# FitAndPheno_plot, ROI, VI, and Rdata
folders <- c("FitAndPheno_plots","ROI","VI", "Rdata") #"./CSV",
for(i in 1:length(folders)){
  if(folders[i] %in% list.dirs("Output",full.names = FALSE)){ 
    #dont do anything if that folder exists in the working directory
  }else{
    #but if its not in the working directory, create it
    dir.create(paste0(OutputDir,folders[i]))
  }
}

#create a site-specific folder in the ROI folder.
#if it already exists, a warning message comes up and nothing happens (nothing overwritten, etc)
dir.create(paste0(OutputDir,"ROI/",Site))

###--- Loop through to draw ROI per ref image #####
# draw with locator. press escape after you finish picking the vertices per image.be patient after pressing escape.
for(i in 1:length(useForROI)){
  site_nameLength<-nchar(strsplit(useForROI[i],"_")[[1]][1])
  roi.name <- substr(useForROI[i],1,site_nameLength) # must be length of nroi; 

  #print(roi.name)
  print(useForROI[i])
  
  #determine if ROI has already been done. if not,
  ROIsdrawn<-DrawMULTIROI(paste0(locOfImgs,Site,"/",Dplymnt[i],"/",useForROI[i]), # the ref picture to use for ROI
          paste0(OutputDir,"ROI/",Site), #path for where to store the copied img with ROI and ROI coords
          nroi = nroi,  file.type = ".JPG",roi.names = paste0(roi.name,"_",Dplymnt[i],"_roi",seq(1,nroi,1)))
  
  #get the vertices from each ROI, for checking consistentcy of placement over time
  if(ROIcheck=="Y"){
  #get the vertices from the ROIs and put them in a for 
  vertices<-data.frame(X1=double(),X2=double(),X3=double(),X4=double(),Y1=double(),Y2=double(),Y3=double(),Y4=double())
  
  for(j in 1:nroi){
    k<-nrow(ROIvertices)
    ROIvertices[k+1,3:6]<-ROIsdrawn[[j]]$vertices$x
    ROIvertices[k+1,7:10]<-ROIsdrawn[[j]]$vertices$y
    ROIvertices$ROI[k+1]<-j
    ROIvertices$Dplymnt[k+1]<-Dplymnt[i]
  }
  }
  
  #The generically names ROIroi.data is written in the outputDir; so rename it and move it to the ROI folder
  file.rename(from=paste0(OutputDir,"ROI/",roi.name,"roi.data.RData"),to=paste0(OutputDir,"ROI/",roi.name,"_",nroi,"ROI","_",Dplymnt[i],".RData"))
  file.copy(paste0(OutputDir,"ROI/",roi.name,"_",nroi,"ROI","_",Dplymnt[i],".RData"), paste0(OutputDir,"ROI/",Site),overwrite = TRUE)
  file.remove(paste0(OutputDir,"ROI/",roi.name,"_",nroi,"ROI","_",Dplymnt[i],".RData"))
}  


#if you had said you wanted to check the consistentcy of ROI placement over time,
if(ROIcheck=="Y"){
  #plot the ROIs across Dplymnts in a single image
  png(filename = paste0(OutputDir,"/ROI/",Site,"/",Site,"_ROIs_axDplymnts",".png"))
  DplymntColors<-rainbow(length(Dplymnt))
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  plot(1, type="n", main=paste0(Site," ROIs over Dplymnts"),
       xlim=c(0, 1), xlab="y",
       ylim=c(0, 0.1+max(ROIvertices[,7:10])), ylab="")
  for(l in 1:nrow(ROIvertices)){
  polygon(x=ROIvertices[l,3:6], y=ROIvertices[l,7:10],border = DplymntColors[which(newDplymntLevels%in%ROIvertices$Dplymnt[l])],lty=ROIvertices$ROI[l])
  }
  
  legend("topright",inset=c(-0.25,0), title="ROI", legend=seq(1,nroi,1),lty=seq(1,nroi,1), 
         col=c("black"), horiz=FALSE,bty="n")
  
  legend("right",inset=c(-0.3,-.3), title="Dplymnt",legend=newDplymntLevels,lty=rep(1,length(Dplymnt)),
         col=DplymntColors, horiz=FALSE,bty="n",cex=0.7)
  
  dev.off()
}
