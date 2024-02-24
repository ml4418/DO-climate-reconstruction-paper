########################################################################################
########################################################################################
#######################                                           ######################
#######################       1. ACER climate reconstructions     ######################
#######################                                           ######################
########################################################################################
########################################################################################

rm(list=ls())

wd<-"D:/PhD Project/DO climate paper/Data and codes/"

if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}
source("D:/PhD Project/DO climate paper/Data and codes/DO climate paper_functions.R", encoding = 'UTF-8')


##############################################################################
####################### Pre-process the modern pollen record
##############################################################################

#Import dataset
smpdsv2_pollen_counts_amalgamated <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_pollen_counts_amalgamated.csv",sep=""),encoding = "UTF-8")
smpdsv2_metadata <- read.csv(paste(wd,"Input data/Pollen/Modern pollen/smpdsv2_metadata.csv",sep=""),encoding = "UTF-8")
colnames(smpdsv2_pollen_counts_amalgamated)[1]<-"ID_SAMPLE"
colnames(smpdsv2_metadata)[1]<-"ID_SITE"

#Get information
unique<-unique(smpdsv2_metadata[,c("longitude","latitude","elevation")])
nrow(smpdsv2_metadata);nrow(unique)

#Merge climate data and pollen data by ID_sample
modern_data<-merge(x=smpdsv2_metadata,y=smpdsv2_pollen_counts_amalgamated,by="ID_SAMPLE",by.y="ID_SAMPLE")
modern_data<-na.omit(modern_data,cols = c("mtco", "mtwa","gdd0","mi")) #remove rows with NAs

taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]
modern_taxa<-modern_taxa/rowSums(modern_taxa)

modern_data_aggregate<-cbind.data.frame(modern_data[,c("longitude","latitude","elevation","mtco", "mtwa","gdd0","mi")],modern_taxa)
modern_data_aggregate<-aggregate(data=modern_data_aggregate,.~longitude+latitude+elevation,FUN=mean)#aggregate duplicate rows
write.csv(modern_data_aggregate,paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat elv.csv",sep=""))


##############################################################################
####################### Training
##############################################################################

modern_data<-read.csv(paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat elv.csv",sep=""),row.names=1)
modern_data$alpha<-get_alpha_from_mi(modern_data$mi)
taxaColMin <- which(colnames(modern_data) == "Abatia")
taxaColMax <- which(colnames(modern_data) == "Zygophyllaceae")
modern_taxa <- modern_data[, taxaColMin:taxaColMax]

#combine Quercus
modern_taxa$Quercus.combined<-modern_taxa$Quercus+modern_taxa$Quercus.deciduous+modern_taxa$Quercus.evergreen
modern_taxa[,c("Quercus","Quercus.deciduous","Quercus.evergreen")]<-NULL

#tidy modern data
modern_taxa <- modern_taxa[, which(colSums(modern_taxa>0)>=10)] #remove taxa less than 10 occurrences
modern_taxa <- modern_taxa[rowSums(modern_taxa)>0,] #remove rows with no taxa information
modern_taxa<-modern_taxa/rowSums(modern_taxa) #re-sum to 1
modern_taxa<-modern_taxa[,sort(colnames(modern_taxa))] #re-order column by names

#training
fit_mtco <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$mtco, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_mtwa <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$mtwa, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.02)
fit_alpha <- fxTWAPLS::TWAPLS.w2(modern_taxa, modern_data$alpha, nPLS = 5, usefx = TRUE, fx_method="pspline",bin=0.002)

fxTWAPLS::plot_train(fit_mtco,3);fxTWAPLS::plot_residuals(fit_mtco,3)
fxTWAPLS::plot_train(fit_mtwa,3);fxTWAPLS::plot_residuals(fit_mtwa,3)
fxTWAPLS::plot_train(fit_alpha,3);fxTWAPLS::plot_residuals(fit_alpha,3)

#############################################################################################
############################  Cross validation 
#############################################################################################

#This part takes too long so it's run on HPC
if(!require(foreach)){install.packages("foreach");library(foreach)}
if(!require(doParallel)){install.packages("doParallel");library(doParallel)}
`%>%` <- magrittr::`%>%`

CPUS<-detectCores()-1

point <- modern_data[, c("longitude", "latitude")]
dist <- fxTWAPLS::get_distance(point, cpus = CPUS)%>% fxTWAPLS::pb()
write.csv(dist, paste(wd,"Output data/ACER reconstructions/Cross validation/distance.csv",sep=""))

# Pseudo removed leave out cross validation
dist<-read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/distance.csv",sep=""),row.names=1)
gc()
options(future.globals.maxSize = 3 * 1024^3)

pseudo_mtco<-fxTWAPLS::get_pseudo(dist, modern_data$mtco, cpus = CPUS) %>% fxTWAPLS::pb()  
pseudo_mtwa<-fxTWAPLS::get_pseudo(dist, modern_data$mtwa, cpus = CPUS) %>% fxTWAPLS::pb()  
pseudo_alpha<-fxTWAPLS::get_pseudo(dist, modern_data$alpha, cpus = CPUS) %>% fxTWAPLS::pb()  


#tf2 pspline
if(!require(foreach)){install.packages("foreach");library(foreach)}

cv_mtco <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$mtco,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_mtco,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 0.02,
                             cpus = CPUS,
                             test_mode = FALSE)  %>% fxTWAPLS::pb()  

cv_mtwa <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$mtwa,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_mtwa,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 0.02,
                             cpus = CPUS,
                             test_mode = FALSE)   %>% fxTWAPLS::pb()  

cv_alpha <- fxTWAPLS::cv.pr.w(modern_taxa,
                             modern_data$alpha,
                             nPLS = 10,
                             fxTWAPLS::TWAPLS.w2,
                             fxTWAPLS::TWAPLS.predict.w,
                             pseudo_alpha,
                             usefx = TRUE,
                             fx_method = "pspline",
                             bin = 0.002,
                             cpus = CPUS,
                             test_mode = FALSE)   %>% fxTWAPLS::pb()  

######### Check the cross validation result
#random test
cv_mtco <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_mtco with Quercus combined.csv",sep=""), row.names=1)
cv_mtwa <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_mtwa with Quercus combined.csv",sep=""), row.names=1)
cv_alpha <- read.csv(paste(wd,"Output data/ACER reconstructions/Cross validation/cv_alpha with Quercus combined.csv",sep=""), row.names=1)

rand_mtco<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_mtco))
rand_mtwa<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_mtwa))
rand_alpha<-as.data.frame(fxTWAPLS::rand.t.test.w(cv_alpha))

rand_mtco$name<-"mtco"
rand_mtwa$name<-"mtwa"
rand_alpha$name<-"alpha"

rand<-rbind.data.frame(rand_mtco[1:5,],rand_mtwa[1:5,],rand_alpha[1:5,])
write.csv(rand, paste(wd,"Output data/ACER reconstructions/Cross validation/random t test result for climate with Quercus combined.csv",sep=""))

##############################################################################
####################### Reconstruction using ACER fossil pollen
##############################################################################
nsig_mtco<-3
nsig_mtwa<-3
nsig_alpha<-4

################ Load fossil pollen and get it formatted

setwd(paste(wd,"Input data/Pollen/Fossil pollen/",sep=""))
listcsv <- dir(pattern = "*.csv") # creates the list of all the csv files in the directory

allcore<-data.frame()
for(k in 1:length(listcsv)){
  
  #load
  fossil_pollen<-read.csv(listcsv[k])
  
  #tidy
  each_core0<-fossil_pollen[,11:ncol(fossil_pollen)] #get only taxa information
  each_core0[is.na(each_core0)]<-0 #replace NA with 0
  each_core<-each_core0
  
  #combine Quercus
  col_Quercus<-grep("Quercus",colnames(each_core))
  if(length(col_Quercus)==1){
    each_core$Quercus.combined<-each_core[,col_Quercus]
  }else{
    each_core$Quercus.combined<-rowSums(each_core[,col_Quercus])
  }
  each_core[,col_Quercus]<-NULL
  
  #get it formatted as modern pollen
  taxa_to_delete<-colnames(each_core)[!(colnames(each_core)%in% colnames(modern_taxa))]
  each_core[,taxa_to_delete]<-NULL
  taxa_to_add<-colnames(modern_taxa)[!(colnames(modern_taxa)%in% colnames(each_core))]
  each_core[,taxa_to_add]<-0
  each_core<-each_core[,order(colnames(each_core))]
  each_core<-each_core/rowSums(each_core)
  
  each_core<-cbind.data.frame(fossil_pollen[,1:10],each_core)
  
  allcore<-rbind.data.frame(allcore,each_core)
}
core<-allcore[,11:ncol(allcore)]

############################# Reconstruction

fossil_mtco<-fxTWAPLS::TWAPLS.predict.w(fit_mtco,core)
fossil_mtwa<-fxTWAPLS::TWAPLS.predict.w(fit_mtwa,core)
fossil_alpha<-fxTWAPLS::TWAPLS.predict.w(fit_alpha,core)

#Get the sample specific errors, use nboot=1000
`%>%` <- magrittr::`%>%`
sse_mtco<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                               modern_climate=modern_data$mtco,
                               fossil_taxa=core,
                               trainfun=fxTWAPLS::TWAPLS.w2,
                               predictfun=fxTWAPLS::TWAPLS.predict.w,
                               nboot=1000,
                               nPLS=5,
                               nsig=nsig_mtco,
                               usefx=TRUE,
                               fx_method = "pspline",
                               bin=0.02,
                               cpus = 6) %>% fxTWAPLS::pb()
sse_mtwa<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                               modern_climate=modern_data$mtwa,
                               fossil_taxa=core,
                               trainfun=fxTWAPLS::TWAPLS.w2,
                               predictfun=fxTWAPLS::TWAPLS.predict.w,
                               nboot=1000,
                               nPLS=5,
                               nsig=nsig_mtwa,
                               usefx=TRUE,
                               fx_method = "pspline",
                               bin=0.02,
                               cpus = 6) %>% fxTWAPLS::pb()
sse_alpha<-fxTWAPLS::sse.sample(modern_taxa=modern_taxa,
                                modern_climate=modern_data$alpha,
                                fossil_taxa=core,
                                trainfun=fxTWAPLS::TWAPLS.w2,
                                predictfun=fxTWAPLS::TWAPLS.predict.w,
                                nboot=20,
                                nPLS=5,
                                nsig=nsig_alpha,
                                usefx=TRUE,
                                fx_method = "pspline",
                                bin=0.002,
                                cpus = 6) %>% fxTWAPLS::pb()
#Use the last significant number of components
core_sig<-cbind.data.frame(allcore[,1:10],
                           fossil_mtco[["fit"]][,nsig_mtco],sse_mtco,
                           fossil_mtwa[["fit"]][,nsig_mtwa],sse_mtwa,
                           fossil_alpha[["fit"]][,nsig_alpha],sse_alpha)
colnames(core_sig)[11:ncol(core_sig)]<-c("mtco","sse_mtco","mtwa","sse_mtwa","alpha","sse_alpha")

#Remove rows with wrong ages -9999
core_sig[core_sig==-9999]=NA 
core_sig<-na.omit(core_sig)

write.csv(core_sig,paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_climate.csv",sep=""))

################################ Add in site information
recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/core_sig_climate.csv",sep=""),row.names = 1)
recon$Tmean<-(recon$mtwa+recon$mtco)/2
recon$sse_Tmean<-sqrt(recon$sse_mtwa^2+recon$sse_mtco^2)/2

siteinfo<-read.csv(paste(wd,"Input data/Pollen/fossil site information.csv",sep=""),row.names = 1)
unique(recon$site_name)

#correct the spelling of names
recon[grep("Abric Roman",recon$site_name),"site_name"]<-"Abric Roman"
recon[which(recon$site_name=="Ca\xe7o"),"site_name"]<-"Cao"
recon[which(recon$site_name=="Col\xf4nia"),"site_name"]<-"Colnia"
recon[which(recon$site_name=="F\xfcramoos"),"site_name"]<-"Framoos"
recon[which(recon$site_name=="Navarr\xe9s"),"site_name"]<-"Navarrs"

#Check whether the site names are all in the file
unique(recon$site_name)%in%unique(siteinfo$site_name)

#Import the site information into the reconstruction file
for(i in 1:nrow(recon)){
  tryCatch({
    name<-recon[i,"site_name"]
    recon[i,c("lon","lat","elv")]<-siteinfo[which(siteinfo$site_name==name),c("long","lat","elev")]
  }, error=function(e){})
}

#Check the age
plot((recon$CLAM_max95+recon$CLAM_min95)/2~recon$CLAM_best);abline(a=0,b=1)
recon$age<-recon$CLAM_best

#Extract only 50~30ka
recon<-recon[which(recon$age>=30000&recon$age<=50000),]
write.csv(recon,paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate.csv",sep=""))

#Get resolution
reso<-get_average_resolution(recon[,c("site_id","age")])
round(mean(reso$avg_reso,na.rm=T),digits=0)

#Get site information of recon_climate
keep<-unique(recon[,c("site_name","site_id")])
siteinfo_to_keep<-merge(keep,siteinfo,by="site_name")
sum(siteinfo_to_keep$site_type=="TERR");sum(siteinfo_to_keep$site_type=="MARI")
write.csv(siteinfo_to_keep,paste(wd,"Output data/ACER reconstructions/siteinfo_50to30ka.csv",sep=""))

########################################################################################
########################################################################################
#######################                                           ######################
#######################       2. Correct alpha for CO2 effect     ######################
#######################                                           ######################
########################################################################################
########################################################################################
rm(list=ls())

wd<-"D:/PhD Project/DO climate paper/Data and codes/"

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(codos)){ install.packages("codos");library(codos)}
source("D:/PhD Project/DO climate paper/Data and codes/DO climate paper_functions.R", encoding = 'UTF-8')

################################### Load past CO2 concentrations

#import ice core data
CO2_source <- read.csv(paste(wd,"Input data/CO2/CO2_WAIS.csv",sep=""))

#import the file to convert GICC05 to AICC2012
#AICC2012 is tuned to GICC05 over the last 60 ka
transfer <- read.csv(paste(wd,"Input data/CO2/GICC05toAICC2012.csv",sep=""))

#remove NAs
CO2_AICC<-na.omit(CO2_source)

#convert timescale
CO2_AICC$AICC2012_BP1950<-as.vector(convert(CO2_AICC$age_calBP/1.0063,transfer))
CO2<-CO2_AICC[,c("AICC2012_BP1950", "CO2_blank_gravity_corrected","CO2_se_1sigma")];colnames(CO2)<-c("age","conc","err_conc")

#bin in 25 years
CO2$age<-ceiling(CO2$age/25)*25 
mean_CO2<-aggregate(data=CO2[,c("age","conc")],.~age,FUN=function(x) mean(x), na.action = na.omit)
err_CO2<-aggregate(data=CO2[,c("age","err_conc")],.~age,FUN=function(x) sqrt(sum(x^2))/length(x), na.action = na.omit)
CO2_binned<-merge(x=mean_CO2,y=err_CO2,by="age")

#add fake points
CO2_fake <- fake(age=as.matrix(CO2_binned$age),
                 value=as.matrix(CO2_binned$conc),
                 err_value=as.matrix(CO2_binned$err_conc),
                 agebin=25)
colnames(CO2_fake)<-c("age","CO2_ppm","err_CO2_ppm")

################################## Load modern temperature

modern_tmp <-read.csv(paste(wd,"Input data/CO2/smpdsv2_climate_reconstructions_tmp.csv",sep=""))
site<-unique(modern_tmp[,c("longitude","latitude")])
plot(site$latitude~site$longitude)
tmp_data<-modern_tmp[,paste("T",1:365,sep="")] #get the temperature data only
tmp_data[tmp_data<=0]<-NA #only keep temperature >0 degree to calculate the mean growing season temperature

modern_tmp$modern_mtgr<-rowMeans(tmp_data,na.rm = TRUE)
modern_tmp$lon<-ceiling(modern_tmp$longitude/5.625)*5.625 #bin in 5.625 degree
modern_tmp$lat<-ceiling(modern_tmp$latitude/5.625)*5.625 #bin in 5.625 degree
modern_mtgr<-aggregate(data=modern_tmp[,c("lon","lat","modern_mtgr")],.~lon+lat,FUN=function(x) mean(x)) #get the mean in each bin

################################ Correct alpha

recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate.csv",sep=""),row.names = 1)
recon$lon<-ceiling(recon$lon/5.625)*5.625 #bin in 5.625 degree
recon$lat<-ceiling(recon$lat/5.625)*5.625 #bin in 5.625 degree
recon$age<-ceiling(recon$age/25)*25 #bin in 25 years
mean<-aggregate(data=recon[,c("lon","lat","elv","age","site_id","mtco","mtwa","Tmean","alpha")],.~lon+lat+elv+site_id+age,FUN=function(x) mean(x)) #get the mean in each bin
err<-aggregate(data=recon[,c("lon","lat","elv","age","site_id","sse_mtco","sse_mtwa","sse_Tmean","sse_alpha")],.~lon+lat+elv+site_id+age,FUN=function(x) sqrt(sum(x^2))/length(x))#get the err in each bin
recon<-merge(mean,err,by=c("lon","lat","elv","age","site_id"))

#approximate mean growing season temperature by MTCO and MTWA
for(i in 1:nrow(recon)){
  ts_tmp<-int_sin(minv=recon[i,"mtco"],maxv=recon[i,"mtwa"],period = 12,period_to_interpolate=c(1:12))
  ts_tmp[ts_tmp<=0]<-NA #only keep temperature >0 degree to calculate the mean growing season temperature
  recon[i,"mtgr"]<- mean(ts_tmp,na.rm=T)
}

#merge data
recon<-merge(recon,CO2_fake,by="age")
recon<-merge(recon,modern_mtgr,by=c("lon","lat"))

for(i in 1:nrow(recon)){
  tryCatch({
    
    #value
    recon[i,"mi_found"]<-find_x(recon[i,"alpha"])
    recon[i,"mi_corrected"] <- codos::corrected_mi(Tc0=recon[i,"modern_mtgr"], 
                                                   Tc1=recon[i,"mtgr"], 
                                                   MI=recon[i,"mi_found"], 
                                                   ca0=340, 
                                                   ca1=recon[i,"CO2_ppm"])
    recon[i,"alpha_corrected"]<-get_alpha_from_mi(recon[i,"mi_corrected"])
    
    #upper
    recon[i,"mi_upper_found"]<-find_x(recon[i,"alpha"]+recon[i,"sse_alpha"])
    recon[i,"mi_upper_corrected"] <- codos::corrected_mi(Tc0=recon[i,"modern_mtgr"], 
                                                         Tc1=recon[i,"mtgr"], 
                                                         MI=recon[i,"mi_upper_found"], 
                                                         ca0=340, 
                                                         ca1=recon[i,"CO2_ppm"])
    recon[i,"alpha_upper_corrected"]<-get_alpha_from_mi(recon[i,"mi_upper_corrected"])
    
    #lower
    recon[i,"mi_lower_found"]<-find_x(recon[i,"alpha"]-recon[i,"sse_alpha"])
    recon[i,"mi_lower_corrected"] <- codos::corrected_mi(Tc0=recon[i,"modern_mtgr"], 
                                                         Tc1=recon[i,"mtgr"], 
                                                         MI=recon[i,"mi_lower_found"], 
                                                         ca0=340, 
                                                         ca1=recon[i,"CO2_ppm"])
    recon[i,"alpha_lower_corrected"]<-get_alpha_from_mi(recon[i,"mi_lower_corrected"])
    
    #sse
    recon[i,"sse_alpha_corrected"]<-(recon[i,"alpha_upper_corrected"]-recon[i,"alpha_lower_corrected"])/2
    
    
  },error = function(e){})
}
write.csv(recon,paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate_binned and CO2 corrected.csv",sep=""))

p<-ggplot(data=recon,aes(alpha,alpha_corrected))+theme_bw()+geom_point(size=0.5)+
  xlim(0,1.26)+ylim(0,1.26)+geom_abline(slope=1,intercept=0,col="red")+
  labs(x=expression(alpha~"before correction"),y=expression(alpha~"after correction"))
ggsave(paste(wd,"Output data/alpha before and after CO2 correction.png",sep=""),p,width=5,height=4.5)

####################################################################################################
####################################################################################################
#######################                                                       ######################
#######################       3. Dynamic time warping to adjust age scale     ######################
#######################                                                       ######################
####################################################################################################
####################################################################################################
rm(list=ls())

wd<-"D:/PhD Project/DO climate paper/Data and codes/"

if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(dtw)){ install.packages("dtw");library(dtw)} #package for dynamic time warping
source("D:/PhD Project/DO climate paper/Data and codes/DO climate paper_functions.R", encoding = 'UTF-8')

#####################################################################

#Load official start dates which are used to make plots
DO_timing1 <- as.matrix(read.csv(paste(wd,"Input data/DO_timing1.csv",sep=""), row.names=1, stringsAsFactors = FALSE))

#Load query
recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate_binned and CO2 corrected.csv",sep=""),row.names = 1)

#Load reference
model<-read.csv(paste(wd,"Input data/LOVECLIM_T_binned_land.csv",sep=""),row.names=1)


#Define the normalization function
norm<-function(x){ (x-mean(x))/sd(x)}

#Adjust
recon_adjusted<-data.frame()

breaks_bin=(DO_timing1[-1]+DO_timing1[-length(DO_timing1)])/2 #use the middle of DO dates to divide the bin
#breaks_bin=c(30000,DO_timing1[5:12],50000)

for(i in unique(recon$site_id)){
  
  each<-recon[which(recon$site_id==i),]
  each<-each[order(each$age),]
  
  if(nrow(each)>2){
    tryCatch(
      {   
        each_adjusted<-data.frame()
        each_fake <- fake(age=each$age,
                          value=each[,c("mtco","mtwa","Tmean","alpha","alpha_corrected")],
                          err_value=each[,c("sse_mtco","sse_mtwa","sse_Tmean","sse_alpha","sse_alpha_corrected")],
                          agebin=25)
        
        lon<-unique(each$lon);lat<-unique(each$lat)
        each_model<-model[which(model$lon==lon&model$lat==lat&
                                  model$age>=min(each$age)& model$age<=max(each$age)),]
        
        #Divide
        each_fake$age_zone<-cut(each_fake$age,breaks=breaks_bin)
        each_model$age_zone<-cut(each_model$age,breaks=breaks_bin)
        
        for(t in unique(each_fake$age_zone)){
          sub_each_fake<-each_fake[which(each_fake$age_zone==t),]
          sub_each_model<-each_model[which(each_model$age_zone==t),]
          
          sub_each_fake$Tmean_norm<-norm(sub_each_fake$Tmean)
          sub_each_model$T_norm<-norm(sub_each_model$T)
          
          #Align
          query<-sub_each_fake$Tmean_norm
          reference<-sub_each_model$T_norm
          alignment<-dtw(query,reference,keep=T)
          #plot(alignment,type="twoway")
          
          wq<-warp(alignment,index.reference=FALSE)
          sub_each_adjusted<-sub_each_fake[wq,]
          sub_each_adjusted[,c("lon","lat","age")]<-sub_each_model[,c("lon","lat","age")]
          
          each_adjusted<-rbind.data.frame(each_adjusted,sub_each_adjusted)
          
        }
        
        each_adjusted$site_id<-unique(each$site_id)
        each_adjusted$elv<-unique(each$elv)

        recon_adjusted<-rbind.data.frame(recon_adjusted,each_adjusted)
        
        p1<-ggplot()+theme_bw()+xlim(30000,50000)+
          geom_point(data=each,aes(age,Tmean),size=0.5)+
          geom_line(data=each,aes(age,Tmean),size=0.3)+
          geom_ribbon(data=each,aes(x=age,y=Tmean,ymin=Tmean-sse_Tmean*1.96,ymax=Tmean+sse_Tmean*1.96),alpha=0.2)+
          geom_point(data=each_adjusted,aes(age,Tmean),col="dodgerblue2",size=0.3)+
          geom_line(data=each_adjusted,aes(age,Tmean),col="dodgerblue2",size=0.3)+
          geom_ribbon(data=each_adjusted,aes(x=age,y=Tmean,ymin=Tmean-sse_Tmean*1.96,ymax=Tmean+sse_Tmean*1.96),alpha=0.2,fill="dodgerblue2")+
          geom_point(data=each_model,aes(age,T),col="red",size=0.3)+
          geom_line(data=each_model,aes(age,T),col="red",size=0.3)+
          geom_ribbon(data=each_model,aes(x=age,y=T,ymin=T-sd_T*1.96,ymax=T+sd_T*1.96),alpha=0.2,fill="red")+
          geom_vline(xintercept=DO_timing1[5:12],linetype="dashed")+
          ggtitle(paste("site_id=",i,"; lon=",lon,"; lat=",lat,"; elv=",unique(each$elv) ))+
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
        
        p2<-ggplot()+theme_bw()+xlim(30000,50000)+
          geom_point(data=each,aes(age,mtco),size=0.5)+
          geom_line(data=each,aes(age,mtco),size=0.3)+
          geom_ribbon(data=each,aes(x=age,y=mtco,ymin=mtco-sse_mtco*1.96,ymax=mtco+sse_mtco*1.96),alpha=0.2)+
          geom_point(data=each_adjusted,aes(age,mtco),col="dodgerblue2",size=0.3)+
          geom_line(data=each_adjusted,aes(age,mtco),col="dodgerblue2",size=0.3)+
          geom_ribbon(data=each_adjusted,aes(x=age,y=mtco,ymin=mtco-sse_mtco*1.96,ymax=mtco+sse_mtco*1.96),alpha=0.2,fill="dodgerblue2")+
          geom_vline(xintercept=DO_timing1[5:12],linetype="dashed")+
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
        
        p3<-ggplot()+theme_bw()+xlim(30000,50000)+
          geom_point(data=each,aes(age,mtwa),size=0.5)+
          geom_line(data=each,aes(age,mtwa),size=0.3)+
          geom_ribbon(data=each,aes(x=age,y=mtwa,ymin=mtwa-sse_mtwa*1.96,ymax=mtwa+sse_mtwa*1.96),alpha=0.2)+
          geom_point(data=each_adjusted,aes(age,mtwa),col="dodgerblue2",size=0.3)+
          geom_line(data=each_adjusted,aes(age,mtwa),col="dodgerblue2",size=0.3)+
          geom_ribbon(data=each_adjusted,aes(x=age,y=mtwa,ymin=mtwa-sse_mtwa*1.96,ymax=mtwa+sse_mtwa*1.96),alpha=0.2,fill="dodgerblue2")+
          geom_vline(xintercept=DO_timing1[5:12],linetype="dashed")+
          theme(axis.title.x = element_blank(), axis.text.x = element_blank())
        
        p4<-ggplot()+theme_bw()+xlim(30000,50000)+
          geom_point(data=each,aes(age,alpha_corrected),size=0.5)+
          geom_line(data=each,aes(age,alpha_corrected),size=0.3)+
          geom_ribbon(data=each,aes(x=age,y=alpha_corrected,ymin=alpha_corrected-sse_alpha_corrected*1.96,ymax=alpha_corrected+sse_alpha_corrected*1.96),alpha=0.2)+
          geom_point(data=each_adjusted,aes(age,alpha_corrected),col="dodgerblue2",size=0.3)+
          geom_line(data=each_adjusted,aes(age,alpha_corrected),col="dodgerblue2",size=0.3)+
          geom_ribbon(data=each_adjusted,aes(x=age,y=alpha_corrected,ymin=alpha_corrected-sse_alpha_corrected*1.96,ymax=alpha_corrected+sse_alpha_corrected*1.96),alpha=0.2,fill="dodgerblue2")+
          geom_vline(xintercept=DO_timing1[5:12],linetype="dashed")
        
        p<-ggarrange(p1,p2,p3,p4,ncol=1)
        
        #The black line is the reconstruction with original age scale, the red line is the LOVECLIM simulation, the blue line is the reconstruction with adjusted age scale 
        ggsave(paste(wd,"Output data/ACER reconstructions/Reconstruction figure/site_id=",i,".png",sep=""),p,width=8,height=12)
        
      },error = function(e){})
  }
  
}

rownames(recon_adjusted)<-1:nrow(recon_adjusted)
recon_adjusted<-na.omit(recon_adjusted)
write.csv(recon_adjusted,paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate_binned and CO2 corrected_age adjusted.csv",sep=""))


####################################################################################################
####################################################################################################
#######################                                                       ######################
#######################             4. Analysis, figures and tables           ######################
#######################                                                       ######################
####################################################################################################
####################################################################################################
rm(list=ls())

wd<-"D:/PhD Project/DO climate paper/Data and codes/"
if(!require(ggplot2)){ install.packages("ggplot2");library(ggplot2)}
if(!require(ggmap)){ install.packages("ggmap");library(ggmap)}
if(!require(ggsn)){ install.packages("ggsn");library(ggsn)}
if(!require(maps)){ install.packages("maps");library(maps)}
if(!require(mapdata)){ install.packages("mapdata");library(mapdata)}
if(!require(egg)){ install.packages("egg");library(egg)}
if(!require(fxTWAPLS)){ install.packages("fxTWAPLS");library(fxTWAPLS)}

world <- map_data("world") 

source("D:/PhD Project/DO climate paper/Data and codes/DO climate paper_functions.R", encoding = 'UTF-8')


################### SMPDSv2 and ACER sites figure

modern_data<-read.csv(paste(wd,"Input data/Pollen/Modern pollen/modern_pollen_aggregated by lon lat elv.csv",sep=""),row.names=1)
recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate.csv",sep=""),row.names = 1)

unique_modern<-unique(modern_data[,c("longitude","latitude","elevation")]);nrow(unique_modern)
unique_fossil<-unique(recon$site_id);length(unique_fossil)

xlab<-c(expression("90"~degree~W),expression("0"~degree~E),expression("90"~degree~E))
ylab<-c(expression("60"~degree~S),expression("0"~degree~N),expression("60"~degree~N))

p1<-ggplot()+theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='black') +
  scale_x_continuous(breaks = c(-90,0,90),labels=xlab)+
  scale_y_continuous(breaks = c(-60,0,60),labels=ylab)+
  geom_point(data=recon,aes(lon,lat),color="red")+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank(),
        axis.text.x=element_blank())

p2<-ggplot()+theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0,color='black') +
  scale_x_continuous(breaks = c(-90,0,90),labels=xlab)+
  scale_y_continuous(breaks = c(-60,0,60),labels=ylab)+
  geom_point(data=modern_data,aes(longitude,latitude),color="red",size=0.5)+
  theme(axis.title.x=element_blank(),axis.title.y=element_blank())

p<-ggarrange(p1,p2,labels=c("(a)","(b)"))

ggsave(paste(wd,"Output data/ACER and SMPDSv2 sites.png",sep=""),p,width=6,height=7)


################################################ changes during D-O events

recon<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate_binned and CO2 corrected.csv",sep=""),row.names = 1)
recon_adjusted<-read.csv(paste(wd,"Output data/ACER reconstructions/Reconstruction/recon_climate_binned and CO2 corrected_age adjusted.csv",sep=""),row.names = 1)
DO_timing1 <- as.matrix(read.csv(paste(wd,"Input data/DO_timing1.csv",sep=""), row.names=1, stringsAsFactors = FALSE))

all_map<-data.frame()

DO_timing1<-ceiling(DO_timing1/25)*25 #bin in 25 years

for(k in 5:12){
  
  t=DO_timing1[k]
  map_win<-recon_adjusted[which(recon_adjusted$age>=(t-600)&recon_adjusted$age<=(t+300)),]
  
  map_k<-data.frame()
  
  for(i in unique(map_win$site_id)){
    sub_map_win<-map_win[map_win$site_id==i,]
    
      lm_mtco<-lm(sub_map_win$mtco ~ poly(sub_map_win$age, 3, raw=TRUE));sub_map_win$mtco_fitted<-lm_mtco[["fitted.values"]]
      lm_mtwa<-lm(sub_map_win$mtwa ~ poly(sub_map_win$age, 3, raw=TRUE));sub_map_win$mtwa_fitted<-lm_mtwa[["fitted.values"]]
      lm_alpha_corrected<-lm(sub_map_win$alpha_corrected ~ poly(sub_map_win$age, 3, raw=TRUE));sub_map_win$alpha_corrected_fitted<-lm_alpha_corrected[["fitted.values"]]
      
      #p1<-ggplot()+theme_bw()+
       #geom_point(data=sub_map_win,aes(age,mtco))+geom_line(data=sub_map_win,aes(age,mtco))+
       #geom_point(data=sub_map_win,aes(age,mtco_fitted),col="red")+geom_line(data=sub_map_win,aes(age,mtco_fitted),col="red")+
       #theme(axis.text.x = element_blank(),axis.title.x = element_blank())+geom_vline(xintercept=t,linetype="dashed")
      #p2<-ggplot()+theme_bw()+
       #geom_point(data=sub_map_win,aes(age,mtwa))+geom_line(data=sub_map_win,aes(age,mtwa))+
       #geom_point(data=sub_map_win,aes(age,mtwa_fitted),col="red")+geom_line(data=sub_map_win,aes(age,mtwa_fitted),col="red")+
       #theme(axis.text.x = element_blank(),axis.title.x = element_blank())+geom_vline(xintercept=t,linetype="dashed")
      #p3<-ggplot()+theme_bw()+
       #geom_point(data=sub_map_win,aes(age,alpha_corrected))+geom_line(data=sub_map_win,aes(age,alpha_corrected))+
       #geom_point(data=sub_map_win,aes(age,alpha_corrected_fitted),col="red")+geom_line(data=sub_map_win,aes(age,alpha_corrected_fitted),col="red")+
       #geom_vline(xintercept=t,linetype="dashed")
      #p<-ggarrange(p1,p2,p3)
      
      change_each_site_DO<-cbind.data.frame(k=k,
                                            site_id=i,
                                            lon=unique(sub_map_win$lon),
                                            lat=unique(sub_map_win$lat),
                                            elv=unique(sub_map_win$elv),
                                            mtco_change=find_change(sub_map_win,lm_mtco,"mtco",t)$change,
                                            mtwa_change=find_change(sub_map_win,lm_mtwa,"mtwa",t)$change,
                                            alpha_corrected_change=find_change(sub_map_win,lm_alpha_corrected,"alpha_corrected",t)$change,
                                            err_mtco_change=find_change(sub_map_win,lm_mtco,"mtco",t)$err_change,
                                            err_mtwa_change=find_change(sub_map_win,lm_mtwa,"mtwa",t)$err_change,
                                            err_alpha_corrected_change=find_change(sub_map_win,lm_alpha_corrected,"alpha_corrected",t)$err_change)
      map_k<-rbind.data.frame(map_k,change_each_site_DO)
    
    
  }

  #get the resolution between start and end
  sample<-recon[which(recon$age>=(t-600)&recon$age<=(t+300)),]
  count_k<-as.data.frame(table(sample$site_id))
  colnames(count_k)<-c("site_id","count_sample")
  map_k<-merge(map_k,count_k,by="site_id",all.x=TRUE)
  map_k[is.na(map_k$count_sample),"count_sample"]<-0
  
  all_map<-rbind.data.frame(all_map,map_k)

}

all_map$lat_zone<-cut(all_map$lat,breaks = c(-90, -23.5, 23.5, 90))
levels(all_map$lat_zone)<-c("SET","TROP","NET")

#get no change sites
summary(all_map$count_sample)
no_change<-all_map[which(all_map$mtco_change==0&all_map$mtwa_change==0&all_map$alpha_corrected_change==0),]
no_change_low_reso<-no_change[which(no_change$count_sample<=3),]

n_due<-as.data.frame(table(all_map$site_id));colnames(n_due)<-c("site_id","n_due")
n_miss<-as.data.frame(table(no_change$site_id));colnames(n_miss)<-c("site_id","n_miss")
n_miss_low_reso<-as.data.frame(table(no_change_low_reso$site_id));colnames(n_miss_low_reso)<-c("site_id","n_miss_low_reso")

count_DO<-merge(n_due,n_miss,by="site_id",all.x=TRUE)
count_DO<-merge(count_DO,n_miss_low_reso,by="site_id",all.x=TRUE)

count_DO[is.na(count_DO)]<-0
sum(count_DO$n_due);sum(count_DO$n_miss);sum(count_DO$n_miss_low_reso)

#combine with site information
siteinfo<-read.csv(paste(wd,"Output data/ACER reconstructions/siteinfo_50to30ka.csv",sep=""),row.names=1)
siteinfo_count<-merge(siteinfo,count_DO,all.x=TRUE,by="site_id")
siteinfo_count<-siteinfo_count[order(siteinfo_count$site_name),]

write.csv(siteinfo_count,paste(wd,"Output data/ACER reconstructions/siteinfo_DO_count.csv",sep=""))

############################## Plot change maps
library(RColorBrewer)
all_map<-all_map[-which(all_map$mtco_change==0&all_map$mtwa_change==0&all_map$alpha_corrected_change==0),]

#mtco
summary(all_map[,"mtco_change"])
breaks_mtco = c(-10,-5,-2.5,-0.5, 0.5,2.5,5,10)
all_map$mtco_valuefactor <- cut(all_map$mtco_change, breaks=breaks_mtco)
levels(all_map$mtco_valuefactor)<-c("<=-5","-5~-2.5","-2.5~-0.5","-0.5~0.5","0.5~2.5","2.5~5",">5")
cols_mtco=rev(brewer.pal(n = 9, name = "RdBu"))[c(1:3,5,7:9)]
names(cols_mtco)<-levels(all_map$mtco_valuefactor)
data_mtco=all_map[,c("k","lon","lat","elv","mtco_valuefactor","count_sample")]
colnames(data_mtco)<-c("k","lon","lat","elv","valuefactor","count_sample")
p_mtco<-plot_change_map(data=data_mtco,cols=cols_mtco,levels=levels(all_map$mtco_valuefactor),legend_name=expression(Delta~"MTCO"~(degree~C)~"  "))
ggsave(file=paste(wd,"Output data/mtco change map.png",sep=""),p_mtco,width=9,height=12)

#mtwa
summary(all_map[,"mtwa_change"])
breaks_mtwa = c(-12,-5,-2.5,-0.5, 0.5,2.5,5,12)
all_map$mtwa_valuefactor <- cut(all_map$mtwa_change, breaks=breaks_mtwa)
levels(all_map$mtwa_valuefactor)<-c("<=-5","-5~-2.5","-2.5~-0.5","-0.5~0.5","0.5~2.5","2.5~5",">5")
cols_mtwa=rev(brewer.pal(n = 9, name = "RdBu"))[c(1:3,5,7:9)]
names(cols_mtwa)<-levels(all_map$mtwa_valuefactor)
data_mtwa=all_map[,c("k","lon","lat","elv","mtwa_valuefactor","count_sample")]
colnames(data_mtwa)<-c("k","lon","lat","elv","valuefactor","count_sample")
p_mtwa<-plot_change_map(data=data_mtwa,cols=cols_mtwa,levels=levels(all_map$mtwa_valuefactor),legend_name=expression(Delta~"MTWA"~(degree~C)~"  "))
ggsave(file=paste(wd,"Output data/mtwa change map.png",sep=""),p_mtwa,width=9,height=12)

#alpha_corrected
summary(all_map[,"alpha_corrected_change"])
breaks_alpha_corrected = c(-0.32,-0.15,-0.05,-0.01,0.01,0.05,0.15,0.32)
all_map$alpha_corrected_valuefactor <- cut(all_map$alpha_corrected_change, breaks=breaks_alpha_corrected)
levels(all_map$alpha_corrected_valuefactor)<-c("<=-0.15", "-0.15~-0.05", "-0.05~-0.01","-0.01~0.01","0.01~0.05","0.05~0.15",">0.15")
cols_alpha_corrected=brewer.pal(n = 9, name = "BrBG")[c(1:3,5,7:9)]
names(cols_alpha_corrected)<-levels(all_map$alpha_corrected_valuefactor)
data_alpha_corrected=all_map[,c("k","lon","lat","elv","alpha_corrected_valuefactor","count_sample")]
colnames(data_alpha_corrected)<-c("k","lon","lat","elv","valuefactor","count_sample")
p_alpha_corrected<-plot_change_map(data=data_alpha_corrected,cols=cols_alpha_corrected,levels=levels(all_map$alpha_corrected_valuefactor),legend_name=expression(Delta~alpha~"    "))
ggsave(file=paste(wd,"Output data/alpha_corrected change map.png",sep=""),p_alpha_corrected,width=9,height=12)

####################################### plot median changes of all DO events

#### mtco
median_map_mtco<-get_median_change(all_map[,c("k","lon","lat","mtco_change")])
summary(median_map_mtco[,"change"])
breaks_median_mtco = c(-6,-2.5,-0.5, 0.5,2.5,6)

median_map_mtco$valuefactor <- cut(median_map_mtco$change, breaks=breaks_median_mtco)
levels(median_map_mtco$valuefactor)<-c("<=-2.5","-2.5~-0.5","-0.5~0.5","0.5~2.5",">2.5")
cols_median_mtco=rev(brewer.pal(n = 9, name = "RdBu"))[c(1:2,5,8:9)]
names(cols_median_mtco)<-levels(median_map_mtco$valuefactor)

p_median_mtco<-ggplot() + theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0.3,size = 0.001,color='grey30')+
  geom_point(data = median_map_mtco, aes(x = lon, y = lat, fill = valuefactor),shape=21, size=3,position="identity")+
  scale_fill_manual(values=cols_median_mtco,limits=levels(median_map_mtco$valuefactor))+
  labs(fill=expression(Delta~"MTCO"~(degree~C)))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "right",
        legend.title=element_text(size=17),
        legend.text=element_text(size=17))

#### mtwa
median_map_mtwa<-get_median_change(all_map[,c("k","lon","lat","mtwa_change")])
summary(median_map_mtwa[,"change"])
breaks_median_mtwa = c(-2.5,-0.5,-0.25,0.25,0.5,2.5)

median_map_mtwa$valuefactor <- cut(median_map_mtwa$change, breaks=breaks_median_mtwa)
levels(median_map_mtwa$valuefactor)<-c("<=-0.5","-0.5~-0.25","-0.25~0.25","0.25~0.5",">0.5")
cols_median_mtwa=rev(brewer.pal(n = 9, name = "RdBu"))[c(1:2,5,8:9)]
names(cols_median_mtwa)<-levels(median_map_mtwa$valuefactor)

p_median_mtwa<-ggplot() + theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0.3,size = 0.001,color='grey30')+
  geom_point(data = median_map_mtwa, aes(x = lon, y = lat, fill = valuefactor),shape=21, size=3,position="identity")+
  scale_fill_manual(values=cols_median_mtwa,limits=levels(median_map_mtwa$valuefactor))+
  labs(fill=expression(Delta~"MTWA"~(degree~C)))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "right",
        legend.title=element_text(size=17),
        legend.text=element_text(size=17))

#### alpha_corrected
median_map_alpha_corrected<-get_median_change(all_map[,c("k","lon","lat","alpha_corrected_change")])
summary(median_map_alpha_corrected[,"change"])
breaks_median_alpha_corrected = c(-0.18,-0.05,-0.01,0.01,0.05,0.18)

median_map_alpha_corrected$valuefactor <- cut(median_map_alpha_corrected$change, breaks=breaks_median_alpha_corrected)
levels(median_map_alpha_corrected$valuefactor)<-c("<=-0.05","-0.05~-0.01","-0.01~0.01","0.01~0.05",">0.05")
cols_median_alpha_corrected=brewer.pal(n = 9, name = "BrBG")[c(1:2,5,8:9)]
names(cols_median_alpha_corrected)<-levels(median_map_alpha_corrected$valuefactor)

p_median_alpha_corrected<-ggplot() + theme_bw()+
  geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0.3,size = 0.001,color='grey30')+
  geom_point(data = median_map_alpha_corrected, aes(x = lon, y = lat, fill = valuefactor),shape=21, size=3,position="identity")+
  scale_fill_manual(values=cols_median_alpha_corrected,limits=levels(median_map_alpha_corrected$valuefactor))+
  labs(fill=expression(Delta~alpha))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(0,0,0,0.7,"cm"))+
  theme(legend.position = "right",
        legend.title=element_text(size=17),
        legend.text=element_text(size=17))


## Put them together
p<-ggarrange(p_median_mtco,p_median_mtwa,p_median_alpha_corrected,ncol=1)

ggsave(paste(wd,"Output data/median change of all DO events.png",sep=""),
       p,width=8,height=10)

######################### relationship between changes, only use high resolution samples

cols=c("darkgoldenrod3","red","dodgerblue2")

p<-ggplot(data=all_map,aes(mtwa_change,mtco_change,col=lat_zone))+theme_bw()+geom_point(size=1)+
  labs(x=expression(Delta~"MTWA"~(degree~C)),y=expression(Delta~"MTCO"~(degree~C)))+
  scale_color_manual(values=cols)

ggsave(paste(wd,"Output data/change_mtco and change_mtwa.png",sep=""),p,width=6,height=4)

p<-ggplot(data=all_map,aes(mtwa_change,alpha_corrected_change,col=lat_zone))+theme_bw()+geom_point(size=1)+
  labs(x=expression(Delta~"MTWA"~(degree~C)),y=expression(Delta~alpha))+
  scale_color_manual(values=cols)

ggsave(paste(wd,"Output data/change_alpha_corrected and change_mtwa.png",sep=""),p,width=6,height=4)

############################# Deming regression of the relationship

if(!require(deming)){ install.packages("deming");library(deming)}

#MTCO and MTWA
fit_SET<-deming(data=all_map[which(all_map$lat_zone=="SET"),],mtco_change~mtwa_change-1,ystd=err_mtco_change,xstd=err_mtwa_change)
fit_TROP<-deming(data=all_map[which(all_map$lat_zone=="TROP"),],mtco_change~mtwa_change-1,ystd=err_mtco_change,xstd=err_mtwa_change)
fit_NET<-deming(data=all_map[which(all_map$lat_zone=="NET"),],mtco_change~mtwa_change-1,ystd=err_mtco_change,xstd=err_mtwa_change)

b_SET<-fit_SET[["coefficients"]][2];ci95_b_SET<-fit_SET[["ci"]][2,];err_b_SET<-(max(ci95_b_SET)-min(ci95_b_SET))/(2*1.96)
b_TROP<-fit_TROP[["coefficients"]][2];ci95_b_TROP<-fit_TROP[["ci"]][2,];err_b_TROP<-(max(ci95_b_TROP)-min(ci95_b_TROP))/(2*1.96)
b_NET<-fit_NET[["coefficients"]][2];ci95_b_NET<-fit_NET[["ci"]][2,];err_b_NET<-(max(ci95_b_NET)-min(ci95_b_NET))/(2*1.96)

b_MTCO_MTWA<-cbind.data.frame(b=c(b_SET,b_TROP,b_NET),
                              err_b=c(err_b_SET,err_b_TROP,err_b_NET),
                              lower95=c(ci95_b_SET["lower 0.95"],ci95_b_TROP["lower 0.95"],ci95_b_NET["lower 0.95"]),
                              upper95=c(ci95_b_SET["upper 0.95"],ci95_b_TROP["upper 0.95"],ci95_b_NET["upper 0.95"]))
b_MTCO_MTWA<-round(b_MTCO_MTWA,digits=3)
rownames(b_MTCO_MTWA)<-c("SET","TROP","NET")
b_MTCO_MTWA<-b_MTCO_MTWA[c("NET","TROP","SET"),]
write.csv(b_MTCO_MTWA,paste(wd,"Output data/b_MTCO_MTWA.csv",sep=""))

#alpha_corrected and MTWA
fit_SET<-deming(data=all_map[which(all_map$lat_zone=="SET"),],alpha_corrected_change~mtwa_change-1,ystd=err_alpha_corrected_change,xstd=err_mtwa_change)
fit_TROP<-deming(data=all_map[which(all_map$lat_zone=="TROP"),],alpha_corrected_change~mtwa_change-1,ystd=err_alpha_corrected_change,xstd=err_mtwa_change)
fit_NET<-deming(data=all_map[which(all_map$lat_zone=="NET"),],alpha_corrected_change~mtwa_change-1,ystd=err_alpha_corrected_change,xstd=err_mtwa_change)

b_SET<-fit_SET[["coefficients"]][2];ci95_b_SET<-fit_SET[["ci"]][2,];err_b_SET<-(max(ci95_b_SET)-min(ci95_b_SET))/(2*1.96)
b_TROP<-fit_TROP[["coefficients"]][2];ci95_b_TROP<-fit_TROP[["ci"]][2,];err_b_TROP<-(max(ci95_b_TROP)-min(ci95_b_TROP))/(2*1.96)
b_NET<-fit_NET[["coefficients"]][2];ci95_b_NET<-fit_NET[["ci"]][2,];err_b_NET<-(max(ci95_b_NET)-min(ci95_b_NET))/(2*1.96)

b_alpha_corrected_MTWA<-cbind.data.frame(b=c(b_SET,b_TROP,b_NET),
                                         err_b=c(err_b_SET,err_b_TROP,err_b_NET),
                                         lower95=c(ci95_b_SET["lower 0.95"],ci95_b_TROP["lower 0.95"],ci95_b_NET["lower 0.95"]),
                                         upper95=c(ci95_b_SET["upper 0.95"],ci95_b_TROP["upper 0.95"],ci95_b_NET["upper 0.95"]))
b_alpha_corrected_MTWA<-round(b_alpha_corrected_MTWA,digits=3)
rownames(b_alpha_corrected_MTWA)<-c("SET","TROP","NET")
b_alpha_corrected_MTWA<-b_alpha_corrected_MTWA[c("NET","TROP","SET"),]
write.csv(b_alpha_corrected_MTWA,paste(wd,"Output data/b_alpha_corrected_MTWA.csv",sep=""))

