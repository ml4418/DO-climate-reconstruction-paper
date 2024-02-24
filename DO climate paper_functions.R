#Define a function to get alpha from MI
get_alpha_from_mi <- function(x) {
  w<-3
  fai<-1/x
  F<-1+fai-(1+fai^w)^(1/w)
  alpha<-1.26*x*F
  return(alpha)
}

#Define a function to find x for a given y
find_x <- function(y) {
  result <- uniroot(function(x) get_alpha_from_mi(x) - y, lower=10^(-6),upper=10)
  return(result$root)
}

#Define a function to do sinusoidal interpolation based on two values, minv and maxv, the lower and upper bounds/peaks of the function.
#The curve is y=a sin(bx-pi/2)+c, using the period from 0 to 2pi, minv=-a+c,maxv=a+c, therefore, a=(maxv-minv)/2,c=(maxv+minv)/2
int_sin <- function(minv,
                    maxv,
                    period,
                    period_to_interpolate) {
  amplitud <- (maxv - minv) / 2
  mid_point <- (maxv + minv) / 2
  y <- amplitud * sin(2 * pi / period * period_to_interpolate - pi / 2) + mid_point
  return(y)
}

#Define a function to convert GICC05_BP2000 to AICC2012_BP1950
convert<-function(old,transfer){
  old<-as.matrix(old)
  new<-matrix(nrow=nrow(old),ncol=1)
  for(i in 1:nrow(old)){
    if(old[i]<60000){
      new[i]<-old[i]
    }else{
      age_to_adjust<-old[i]
      row_closest<-which.min(abs(transfer$GICC05AgeBP1950 - age_to_adjust))
      age_closest<-as.numeric(transfer[row_closest,"GICC05AgeBP1950"])
      if(age_to_adjust >= age_closest){
        y1<-as.numeric(transfer[row_closest,"AICC2012AgeBP1950"])
        y2<-as.numeric(transfer[row_closest+1,"AICC2012AgeBP1950"])
        x1<-as.numeric(transfer[row_closest,"GICC05AgeBP1950"])
        x2<-as.numeric(transfer[row_closest+1,"GICC05AgeBP1950"])
        age_adjusted<-y1+(y2-y1)*(age_to_adjust-x1)/(x2-x1)
      }else if(age_to_adjust < age_closest){
        y1<-as.numeric(transfer[row_closest-1,"AICC2012AgeBP1950"])
        y2<-as.numeric(transfer[row_closest,"AICC2012AgeBP1950"])
        x1<-as.numeric(transfer[row_closest-1,"GICC05AgeBP1950"])
        x2<-as.numeric(transfer[row_closest,"GICC05AgeBP1950"])
        age_adjusted<-y1+(y2-y1)*(age_to_adjust-x1)/(x2-x1)
      }
      new[i]<-age_adjusted
      
    }
    
  }
  return(new)
}

#Define a function to add fake points
fake<-function(age,value,err_value,agebin){
  data<-cbind.data.frame(age,value,err_value)
  data_new<-data.frame(age=seq(round(min(age),digits=0),max(age),by=agebin))
  
  for(row in 1:nrow(data_new)){
    age_new=data_new[row,"age"]
    row_closet_in_data=which.min(abs(age_new-data[,"age"]))
    age_closet_in_data=data[row_closet_in_data,"age"]
    
    if(age_new==age_closet_in_data){
      value_new=value[row_closet_in_data,]
      err_value_new=err_value[row_closet_in_data,]
      
    }else if(age_new>age_closet_in_data){
      y1<-as.numeric(value[row_closet_in_data,])
      y2<-as.numeric(value[row_closet_in_data+1,])
      err_y1<-as.numeric(err_value[row_closet_in_data,])
      err_y2<-as.numeric(err_value[row_closet_in_data+1,])
      x1<-as.numeric(data[row_closet_in_data,"age"])
      x2<-as.numeric(data[row_closet_in_data+1,"age"])
      value_new<-y1+(y2-y1)*(age_new-x1)/(x2-x1) #y1*(1-(age_new-x1)/(x2-x1))+y2*(age_new-x1)/(x2-x1)
      err_value_new<-sqrt( err_y1^2*(1-(age_new-x1)/(x2-x1))^2 + err_y2^2*((age_new-x1)/(x2-x1))^2 )
      
    }else{
      y1<-as.numeric(value[row_closet_in_data-1,])
      y2<-as.numeric(value[row_closet_in_data,])
      err_y1<-as.numeric(err_value[row_closet_in_data-1,])
      err_y2<-as.numeric(err_value[row_closet_in_data,])
      x1<-as.numeric(data[row_closet_in_data-1,"age"])
      x2<-as.numeric(data[row_closet_in_data,"age"])
      value_new<-y1+(y2-y1)*(age_new-x1)/(x2-x1) #y1*(1-(age_new-x1)/(x2-x1))+y2*(age_new-x1)/(x2-x1)
      err_value_new<-sqrt( err_y1^2*(1-(age_new-x1)/(x2-x1))^2 + err_y2^2*((age_new-x1)/(x2-x1))^2 )
      
    }#end of the if function
    
    if(length(value_new)==0){
      data_new[row,colnames(data)]=NA
    }else{
      data_new[row,colnames(data)]=c(age_new,value_new,err_value_new)
    }   
  }#end of the for loop
  
  return(data_new)
}#overall end

#Define a function to find change at a given site and DO
find_change<-function(sub_map_win,lm,var.name,t){
  a<-lm[["coefficients"]]["poly(sub_map_win$age, 3, raw = TRUE)3"]
  b<-lm[["coefficients"]]["poly(sub_map_win$age, 3, raw = TRUE)2"]
  c<-lm[["coefficients"]]["poly(sub_map_win$age, 3, raw = TRUE)1"]
  d<-lm[["coefficients"]]["(Intercept)"]
  f_polynomial<-function(x){ a*x^3+b*x^2+c*x+d }
  min_t_polynomial <- optimize(f_polynomial, c(t-600, t+300) )$minimum
  max_t_polynomial <- optimize(function(x)-f_polynomial(x), c(t-600, t+300) )$minimum
  min<-min(sub_map_win[which(sub_map_win$age>=min((min_t_polynomial-100),sub_map_win$age)&sub_map_win$age<=max((min_t_polynomial+100),sub_map_win$age)),var.name])
  max<-max(sub_map_win[which(sub_map_win$age>=min((max_t_polynomial-100),sub_map_win$age)&sub_map_win$age<=max((max_t_polynomial+100),sub_map_win$age)),var.name])
  
    if(min_t_polynomial<=max_t_polynomial){
      #decrease from max to min
      t_start_data<-sub_map_win[which(sub_map_win[,var.name]==max),]
      t_end_data<-sub_map_win[which(sub_map_win[,var.name]==min),]
      change<-unique(t_end_data[,var.name])-unique(t_start_data[,var.name])
      err_change<-sqrt(unique(t_end_data[,paste("sse_",var.name,sep="")])^2+unique(t_start_data[,paste("sse_",var.name,sep="")])^2)
    }else{
      #increase from min to max
      t_start_data<-sub_map_win[which(sub_map_win[,var.name]==min),]
      t_end_data<-sub_map_win[which(sub_map_win[,var.name]==max),]
      change<-unique(t_end_data[,var.name])-unique(t_start_data[,var.name])
      err_change<-sqrt(unique(t_end_data[,paste("sse_",var.name,sep="")])^2+unique(t_start_data[,paste("sse_",var.name,sep="")])^2)
    }

  output<-cbind.data.frame(change,err_change)
  return(output)
}

#Define a function to get common legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#Define a function to plot maps for the changes
plot_change_map<-function(data,cols,levels,legend_name){
  
  list<-list()
  
  for(k in 5:12){
    
    map_change<-data[which(data$k==k),]
    #map_change$elv_zone<-as.factor(ifelse(map_change$elv>1500,">1500m","<=1500m"))
    #map_change$count_zone<-as.factor(ifelse(map_change$count_sample<=3,"low reso","high reso"))
    #map_change$combined_factor<-interaction(map_change$elv_zone,map_change$count_zone, sep = "  ")
    #shapes=c(15,17,0,2)
    #names(shapes)=c("<=1500m  high reso",">1500m  high reso","<=1500m  low reso", ">1500m  low reso")
    
    
    p_each<-ggplot() + theme_bw()+
      geom_polygon(data = world[world$region!="Antarctica",],aes(x=long, y = lat, group = group),alpha=0.3,size = 0.001,color='grey30')+
      geom_point(data = map_change, aes(x = lon, y = lat, fill = valuefactor),shape=21, size = 2.5, position="identity")+
      scale_fill_manual(values=cols,limits=levels)+
      #scale_shape_manual(values=shapes,limits=c("<=1500m  high reso",">1500m  high reso","<=1500m  low reso", ">1500m  low reso"))+
      labs(fill=legend_name)+
            theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank(),
            plot.margin = margin(0,0,0,0,"cm"))+
      theme(legend.position = "bottom",
            legend.title = element_text(size=13),legend.text = element_text(size=13))+
      annotate("text", y= 90, x =-170,label=paste("DO",k),size=5)
    
    legend<-get_legend(p_each)
    list[[k-4]]<-p_each
    gc()
  }
  
  p1<-list[[1]]+theme(legend.position = "none")
  p2<-list[[2]]+theme(legend.position = "none")
  p3<-list[[3]]+theme(legend.position = "none")
  p4<-list[[4]]+theme(legend.position = "none")
  p5<-list[[5]]+theme(legend.position = "none")
  p6<-list[[6]]+theme(legend.position = "none")
  p7<-list[[7]]+theme(legend.position = "none")
  p8<-list[[8]]+theme(legend.position = "none")
  
  DO_map<-grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,legend,ncol=2,
                       layout_matrix=rbind(cbind(1,2),cbind(1,2),cbind(3,4),cbind(3,4),
                                           cbind(5,6),cbind(5,6),cbind(7,8),cbind(7,8),9))
  return(DO_map)
}

#Get average age resolution for each site
get_average_resolution<-function(data){
  
  colnames(data)<-c("site_id","age")
  resolution<-data.frame(matrix(nrow=0,ncol=2))
  
  for(i in unique(data$site_id)){
    subdata<-data[which(data$site_id==i),]
    subdata<-subdata[order(subdata$age),]
    resolution_subdata<-cbind.data.frame(i,mean(subdata[-1,"age"]-subdata[-nrow(subdata),"age"]))
    resolution<-rbind.data.frame(resolution,resolution_subdata)
  }
  colnames(resolution)<-c("site_id","avg_reso")
  return(resolution)
  
}

#Get median change of DO events in the signal
get_median_change<-function(data){
  colnames(data)<-c("k","lon","lat","change")
  #data$sign<-ifelse(data$change>=0,1,-1)
  #mean_sign<-aggregate(data[,c("k","lon","lat","sign")],.~lon+lat,FUN=function(x) mean(x))
  #mean_sign$consistency<-ifelse(abs(mean_sign$sign)==1,"consistent","not consistent")
  median_change<-aggregate(data[,c("k","lon","lat","change")],.~lon+lat,FUN=function(x) median(x))
  return(median_change)
}




