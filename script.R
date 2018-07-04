colnames(datnew_male)

CRC.polyps.JNF$call
coef_name=c("t_start0","t_stop0","polyp_event",
            "idnew","famid","CRC",
            "firstvisitage_trans",
           "age_prev_trans",
           "probandage_trans",
           "gender" )
datapred=datnew_male[,coef_name]
pred=prediction(fit = CRC.polyps.JNF, data=datapred,t=0,
           window=1:15,event = "Terminal", conditional = F,
           individual = 0, MC.sample = 500)
