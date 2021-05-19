##############################################################
# Timelapse_Phenopix_4_Organizing Results
#bring in all the r data from the extracted phenology
#organize the fitted curves/predicted values per site per year
#caluclate the DHI Components:
#annual total greenness (ie area under curve);
#annual degree of variation (mean/st.dev), ie., coefficient of variation
# annual maximum greenness (versus orginally the minimum from Coops and Radeloff formulations)
##############################################################


#Packages #####
#these arent all necessary; some are remnant from a previous verson of this script that had more code
library(phenopix)
library(jpeg)
library(rasterImage)
library(lubridate)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(zoo)
require(reshape2)
library(stringr)
library(data.table)
library(DescTools)
library(ggpubr)
library(stringr)
library(tidyr)

###---- Define and Set thingsTHIS NEEDS TO BE SPECIFIC TO YOUR PROJECT AND YOUR FOLDER STRUCTURE. ----####

#where the data are
#InputDir <-"Z:/Camera_Trap_Projects/Active Projects/Algar/3. Data/3.4 Working Projects/CatS/Phenology/"
#setwd(InputDir)

 
#
OutputDir <-"./Output/RData/"

#the sites whose phenology data you want to use now that youve extracted it
# in theory, something like this:
SitesToCollate<-paste0("CATH",str_pad(seq(8,46,1), 2, pad = "0"))
#for the purposes of this demo, this is just
SitesToCollate<-"CATH08"


#create an empty dataframe that the  dates from each site will go into
Metrics_Spline<-data.frame()
Metrics_Klosterman<-data.frame()
Metrics_TRS<-data.frame()

FittedPreds_ROI1_Spline<-list()
FittedPreds_ROI1_Klosterman<-list()
FittedPreds_ROI1_Trs<-list()
Metrics_ROI<-list()

###---- Bring in Data #####
for(i in 1:length(SitesToCollate)){
  print(SitesToCollate[i])
  
  #read in the r data for each site
  realI<-i
  realOutputDir<-OutputDir
  load(paste0("Rdata/PhenopixAnalysisResults_",SitesToCollate[i],".RData"))
  
  i<-realI # redefine i because we dont want to use the 'i' value from the Rdata we just read in
  OutputDir<-realOutputDir
  
  #Get All Metrics_Fitted_AllROIs_overTime values for that site 
  #Rename AllMetrics_Fitted_AllROIs_overTime to be site-specific.  
  assign(paste0(SitesToCollate[i],"_Metrics"),AllMetrics_Fitted_AllROIs_overTime)
  Metrics_ROI[[i]]<-AllMetrics_Fitted_AllROIs_overTime
  names(Metrics_ROI)[i]<-SitesToCollate[i]
  
  #put the metric dates into dataframe
  Metrics_Klosterman <- rbind(Metrics_Klosterman,AllMetrics_Fitted_AllROIs_overTime$KlKl[which(AllMetrics_Fitted_AllROIs_overTime$KlKl$ROI==1),])
  Metrics_Spline <- rbind(Metrics_Spline,AllMetrics_Fitted_AllROIs_overTime$SpKl[which(AllMetrics_Fitted_AllROIs_overTime$SpKl$ROI==1),])
  Metrics_TRS <- rbind(Metrics_TRS,AllMetrics_Fitted_AllROIs_overTime$SpTRS[which(AllMetrics_Fitted_AllROIs_overTime$SpTRS$ROI==1),])
  
  #put the predicted values into lists
  temp<-get(paste0("Fitted_preds_",SitesToCollate[i],"_ROI1_allYrs"))
 
  FittedPreds_ROI1_Spline[[i]]<-temp$SpKl
  if(is.data.frame(temp$SpKl)){
    names(FittedPreds_ROI1_Spline)[i]<-SitesToCollate[i]
   
    }
      
  FittedPreds_ROI1_Klosterman[[i]]<-temp$KlKl
  if(is.data.frame(temp$KlKl)){
    names(FittedPreds_ROI1_Klosterman)[i]<-SitesToCollate[i]
  }
  
  FittedPreds_ROI1_Trs[[i]]<-temp$SpTRS
  if(is.data.frame(temp$SpTRS)){
    names(FittedPreds_ROI1_Trs)[i]<-SitesToCollate[i]
    }
  }  

Metrics_Klosterman$Fit<-rep("Kl",nrow(Metrics_Klosterman))
Metrics_Spline$Fit<-rep("Sp",nrow(Metrics_Spline))
Metrics_TRS$Fit<-rep("TRS",nrow(Metrics_TRS))

head(Metrics_Klosterman)
head(Metrics_Spline)
head(Metrics_TRS)
str(FittedPreds_ROI1_Trs)

#in the above loop, things brought in from the Rdata are rewritten because theyre all by the same name
# so instead of removing them each time in that loop, i've left it so it could just be  be written over each time
# and then just delete once at the end , ie., now, to clean up the environment.
# rm(list=setdiff(ls(), c("Metrics_Klosterman",
#                         "Metrics_Spline",
#                         "Metrics_TRS",
#                         "FittedPreds_ROI1_Spline",
#                         "FittedPreds_ROI1_Klosterman",
#                         "FittedPreds_ROI1_Trs",
#                         "Metrics_ROI",
#                         "OutputDir",
#                         "SiteCovs",
#                         "SitesToCollate","CovsToLookAt")))

FittedPreds_ROI1_Spline<-FittedPreds_ROI1_Spline
# convert things from lists to dataframes
FittedPreds_ROI1_Trs_df <- do.call(rbind, FittedPreds_ROI1_Trs)
FittedPreds_ROI1_Trs_df$Fit<-"Spline_TRS"

FittedPreds_ROI1_Spline_df <- do.call(rbind, FittedPreds_ROI1_Spline)
FittedPreds_ROI1_Spline_df$Fit<-"Spline_Klosterman"

FittedPreds_ROI1_Klosterman_df <- do.call(rbind, FittedPreds_ROI1_Klosterman)
FittedPreds_ROI1_Klosterman_df$Fit<-"Klosterman_Klosterman"

FittedPreds_ROI1_allFits_df<-rbind(FittedPreds_ROI1_Trs_df,FittedPreds_ROI1_Spline_df,FittedPreds_ROI1_Klosterman_df)


####---- explore via plots #####

#plot predicted values

#plot all the fitted curves for all sites together
whichPoints<-which(FittedPreds_ROI1_allFits_df$Date>as.Date("2019-06-30"))
plot_allFits<-ggplot(FittedPreds_ROI1_allFits_df[,], aes(x = Date, y = Pred)) + 
  geom_line(aes(color = Fit,group = Site)) +
  #geom_smooth(method = "loess", size = 1)+
  #stat_summary(aes(group=bucket), fun.y=mean, geom="line", colour="black")+
  facet_grid(rows=vars(Fit),switch="y")+
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="black")+
  theme_minimal()



################

#get means per year
meanPhenoDates_TRS<-Metrics_TRS
meanPhenoDates_TRS$Year<-year(meanPhenoDates_TRS$`50%`)

meanPhenoDates_TRS<-meanPhenoDates_TRS %>% #,
  group_by(Event,Year) %>%
  summarize(meanDate = mean(`50%`))

#turn the phenology event column into a factor
meanPhenoDates_TRS$Event <- factor(meanPhenoDates_TRS$Event, levels = c("Greenup", "Maturity", "Senescence", "Dormancy"))

plot_greenPheno_AllSitesAllYrs<-ggplot(FittedPreds_ROI1_allFits_df[FittedPreds_ROI1_allFits_df$Fit=="Spline_TRS",]) + 
  geom_vline(data= meanPhenoDates_TRS[meanPhenoDates_TRS$Event=="Greenup",],show.legend = FALSE,  aes(xintercept=meanDate,color="Date: Green Up" ),size=1)+
  geom_vline(data= meanPhenoDates_TRS[meanPhenoDates_TRS$Event=="Maturity",], show.legend = FALSE,aes(xintercept=meanDate,color="Date: Maturity"),size=1)+
  geom_vline(data= meanPhenoDates_TRS[meanPhenoDates_TRS$Event=="Senescence",], show.legend = FALSE,aes(xintercept=meanDate,color="Date: Senescence"),size=1)+
  geom_vline(data= meanPhenoDates_TRS[meanPhenoDates_TRS$Event=="Dormancy",], show.legend = FALSE,aes(xintercept=meanDate,color="Date: Dormancy"),size=1)+
  geom_line(aes(x = Date, y = Pred,group = Site,color = "Curve: Per Site"),show.legend = TRUE) +
  stat_summary(fun.y=mean, aes(x = Date, y = Pred,group=1, color="Curve: Mean"), geom="line",show.legend = TRUE)+
  scale_color_manual( values=c( `Date: Green Up` ="green",`Date: Maturity`="darkgreen",
                                `Date: Senescence`="brown",`Date: Dormancy`="lightblue",
                                `Curve: Per Site`="grey",`Curve: Mean`="black"),
                      name = "Extracted Dates and Fitted Curves")+
  #guides(size = guide_legend(title='Fitted Curves'))+
  geom_hline(aes(yintercept=-Inf)) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),strip.placement = "outside")+
  labs(title="Vegetation Phenology",
       x ="Year (time)", y = "Fitted Relative Greenness Index")


#turn the df of mean phenology dates per year across sites from long to wide format
meanPhenoDates_TRS_wide <- as.data.frame(t(pivot_wider(meanPhenoDates_TRS,names_from = c(Year),values_from=meanDate)))
colnames(meanPhenoDates_TRS_wide)<-c("Dormancy" ,"Greenup" ,"Maturity" ,"Senescence")
meanPhenoDates_TRS_wide<-meanPhenoDates_TRS_wide[-1,]
meanPhenoDates_TRS_wide$Dormancy<-as.Date(meanPhenoDates_TRS_wide$Dormancy)
meanPhenoDates_TRS_wide$Greenup<-as.Date(meanPhenoDates_TRS_wide$Greenup)
meanPhenoDates_TRS_wide$Maturity<-as.Date(meanPhenoDates_TRS_wide$Maturity)
meanPhenoDates_TRS_wide$Senescence<-as.Date(meanPhenoDates_TRS_wide$Senescence)

#add a Year column
meanPhenoDates_TRS_wide<-cbind(Year=as.numeric(sapply(strsplit(rownames(meanPhenoDates_TRS_wide),"_"), `[`, 1) ),
                                       meanPhenoDates_TRS_wide)
head(meanPhenoDates_TRS_wide)

#plot the curves per site across all years, with the mean growing season in the background
plot_greenPheno_AllSitesAllYrs_v2<-ggplot(FittedPreds_ROI1_allFits_df[FittedPreds_ROI1_allFits_df$Fit=="Spline_TRS",]) + 
     geom_rect(data=meanPhenoDates_TRS_wide, inherit.aes=FALSE,
                             aes(xmin=meanPhenoDates_TRS_wide$Greenup,
                                                   xmax=meanPhenoDates_TRS_wide$Senescence,
                                                   ymin=-Inf,ymax=Inf,fill="Growng Season"), alpha=0.2,show.legend = FALSE)+ #
     geom_line(aes(x = Date, y = Pred,group = Site,color = "Per Site"),show.legend = TRUE) +
     scale_fill_manual( values=c( `Growng Season` ="lightgreen"),name = "")+
      stat_summary(fun.y=mean, aes(x = Date, y = Pred,group=1, color="Mean"), geom="line",show.legend = TRUE)+
      scale_color_manual( values=c( #`Date: Green Up` ="green",#`Date: Maturity`="darkgreen",
                                    #`Date: Senescence`="brown",#`Date: Dormancy`="lightblue",
                                       `Per Site`="grey",`Mean`="black"),
                           name = "Fitted Spline")+
     #guides(size = guide_legend(title='Fitted Curves'))+
     geom_hline(aes(yintercept=-Inf)) + 
     theme_classic()+
     theme(plot.title = element_text(hjust = 0.5),strip.placement = "outside")+
     labs(#title="Vegetation Phenology (2015-2019) \nAcross 73 sites and 5 restoration categories",
            x ="Year (time)", y = "Fitted Relative Greenness Index")


####--- calculate DHI stats per site and put them in a growing dataframe #####
PhenoStats_perSite_perYr<-data.frame()
theseSites<-seq(1,length(FittedPreds_ROI1_Spline),1)[-which(names(FittedPreds_ROI1_Spline)=="")]

#in this demo, theseSites will be empty because of the [] subsetting by a value of 0, so 
if(length(theseSites)==0){
  theseSites<-1
}

for(j in theseSites){ #per site
  
  print(names(FittedPreds_ROI1_Spline)[j])
  
  #get the years available per site
  PhenoStats_thatSite_perYr<-as.data.frame(table(year(FittedPreds_ROI1_Spline[[j]]$Date))) #days per year
  
  #per year
  for(k in 1:nrow(PhenoStats_thatSite_perYr)){ 
    #get the data points
    dataPts_wiThatYr<-which(year(FittedPreds_ROI1_Spline[[j]]$Date)==PhenoStats_thatSite_perYr$Var1[k]) #find the portion of data within that year
    
    #calculate the area under curve
    auc_thatYr<-AUC(dataPts_wiThatYr,FittedPreds_ROI1_Spline[[j]]$Pred[dataPts_wiThatYr],method="trapezoid")
    PhenoStats_thatSite_perYr$AUC[k]<-auc_thatYr
    
    #write 2 lines similar (almost identical) to above to the max predicted/smoothed value per year
    #ie instead of using the AUC function, use max()..
    max_thatYr<-max(dataPts_wiThatYr,FittedPreds_ROI1_Spline[[j]]$Pred[dataPts_wiThatYr])
    PhenoStats_thatSite_perYr$max[k]<-max_thatYr    
    
    #calculate the coefficient of variation
    cv_thatYr<-mean(FittedPreds_ROI1_Spline[[j]]$Pred[dataPts_wiThatYr])/sd(FittedPreds_ROI1_Spline[[j]]$Pred[dataPts_wiThatYr])
    PhenoStats_thatSite_perYr$CV[k]<-cv_thatYr

  }
  PhenoStats_thatSite_perYr$Site<-names(FittedPreds_ROI1_Spline)[j]
   PhenoStats_perSite_perYr<-rbind(PhenoStats_perSite_perYr,PhenoStats_thatSite_perYr)
}

# give more descriptive column names
colnames(PhenoStats_perSite_perYr)<-c("Year","nDays","AUC","Max Greenness Index","Seasonality","Site")

#keep only for those with at least a minimum proportion of the year 
#because these stats arent reliable if theyre only based on a few days, 

minPropOfYr<-0.05*365
PhenoStats_perSite_perYr<-PhenoStats_perSite_perYr[PhenoStats_perSite_perYr$nDays>minPropOfYr,]
head(PhenoStats_perSite_perYr)


