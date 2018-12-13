

################ Run Pseudodata ##############


files2<-c(
  "Calibration_1", "Calibration_2","Calibration_3",
  "Calibration_4","Calibration_5","Calibration_6", "BangBang_1", "BangBang_2",
  "DynStim_1", "DynStim_2", "DynStim_3", "DynStim_4", "DynStim_5", "DynStim_6", "DynStim_7", 
  "DynStim_8", "DynStim_9", "DynStim_10",
  "DynStim_11", "DynStim_12", "DynStim_13", "DynStim_14", "DynStimTooFast", "PI_1", "PI_2", "PI_3",
  "Constant_1", "Constant_2")


files3 <- c("Calibration_4","Calibration_5","Calibration_6", "BangBang_2",
            "DynStim_2, DynStim_3", "DynStim_9", "DynStim_11")

files4 <- c("Calibration_4_Fake2","Calibration_5_Fake2","Calibration_6_Fake2", "BangBang_2_Fake2",
            "DynStim_2_Fake2, DynStim_3_Fake2", "DynStim_9_Fake2", "DynStim_11_Fake2")

for(x in files3){
  pseudoData(x, "fit_DynStim_2", "Fake2")
}







# for(x in files2){
#   pseudoData(x, "fit_Calibration_6", "Fake3")
# }
# 

