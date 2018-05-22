#=========set up============
if(!require(frailtypack,quietly = TRUE)){
  install.packages("frailtypack")
  library(frailtypack)
}
if(!require(gtools,quietly = TRUE)){
  install.packages("gtools")
  library(gtools)
}


#setwd("M:/yuliang/project1")
#datals=read.csv('datals.csv')
datalsnew=read.csv('datalsnew.csv')


#=======Brier score====================
#Brier score for prediction error part1
datasub = datalsnew
predtime = 1
window = seq(1,25,1)
fwindow = predtime + window
matsum1 = matrix(rep(0,25),nrow = 1)
#split all data into 6 folds and each with its complimentary
#the following is a huge loop

for(i in 1:6){
  # if(i==1){
  #   datapred1 = datasub[datasub$famid %in% c(1:40),]
  #   datasub1 = datasub[!datasub$famid %in% c(1:40),]
  # }
  if(i==2){
    datapred1 = datasub[datasub$famid %in% c(41:80),]
    datasub1 = datasub[!datasub$famid %in% c(41:80),]
  }
  if(i==3){
    datapred1 = datasub[datasub$famid %in% c(81:120),]
    datasub1 = datasub[!datasub$famid %in% c(81:120),]
  }
  if(i==4){
    datapred1 = datasub[datasub$famid %in% c(121:160),]
    datasub1 = datasub[!datasub$famid %in% c(121:160),]
  }
  if(i==5){
    datapred1 = datasub[datasub$famid %in% c(161:200),]
    datasub1 = datasub[!datasub$famid %in% c(161:200),]
  }
  if(i==6){
    datapred1 = datasub[datasub$famid %in% c(201:242),]
    datasub1 = datasub[!datasub$famid %in% c(201:242),]
  }
  #Brier score for prediction error part2
  set.seed(111)#??is this necessary??
  CRC.polyps.JF5 = frailtyPenal(Surv(t_start0,t_stop0,polyp_event)~#polyp_event??
                                  cluster(idnew)+
                                  firstvisitage_trans+
                                  age_prev_trans+
                                  cum_adenoma_prev+#cum_adenoma_prev
                                  probandage_trans+
                                  gender+
                                  terminal(CRC),
                                
                                formula.terminalEvent=
                                  ~firstvisitage_trans+
                                  num_adenoma_beforecrc+
                                  probandage_trans+
                                  gender,
                                
                                recurrentAG = T,
                                n.knots = 5,
                                #kappa = c(1, 1e-10),
                                kappa = c(6.47498, 9.74423e-19),
                                maxit = 40,
                                data=datasub1#different from JF
  )
  
  summary(CRC.polyps.JF5)
  #Extract the time-to-CRC
  datapred1$timetocrc = na.replace(datapred1$timetocrc,100)#replace na    with 100? why??
  timetocrc.vec = aggregate(datapred1[,c('timetocrc')],list(datapred1$idnew),mean)
  #Dynamic prediction using the training data
  #ERROR: Attempting to do predictions with a wrong model??
  #Error starts from i=2
  pred_bivariate3 = prediction(CRC.polyps.JF5,
                               data = datapred1,
                               t=predtime,
                               window = window,
                               event = "Terminal",
                               conditional = F)
  predictions_bivariate3 = pred_bivariate3$pred1
  timeto = as.vector(timetocrc.vec[,2])
  timetocrc.mat = replicate(length(window),timeto)
  fwindow.mat = t(replicate(length(timeto),fwindow))
  obs.mat = fwindow.mat > timetocrc.mat
  # Calculate Brier score
  diff.p2 = (predictions_bivariate3-obs.mat)^2
  data1 = as.data.frame(rep=i, diff.p2)
  colnames(matsum1) = c(colnames(data1))
  matsum1 = rbind.data.frame(matsum1,data1)
}



#====joint frailty model for terminal event CRC and recurrent polyps==========

#?? how to choose the covariates?
CRC.polyps.JF = frailtyPenal(Surv(t_start0,t_stop0,polyp_event)~#polyp_event??
                               cluster(idnew)+
                               firstvisitage_trans+
                               age_prev_trans+
                               cum_adenoma_prev+#cum_adenoma_prev
                               probandage_trans+
                               gender+
                               terminal(CRC),
                             
                             formula.terminalEvent=
                               ~firstvisitage_trans+
                               num_adenoma_beforecrc+
                               probandage_trans+
                               gender,
                             
                             recurrentAG = T,
                             n.knots = 5,
                             kappa = c(6.47498, 9.74423e-19),
                             maxit = 40,
                             data=datalsnew)
summary(CRC.polyps.JF)
