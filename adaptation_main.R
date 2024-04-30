
#----------------------------------------------------------------------


fullDataFrame = read.csv( file="motor_adaptation_archerfish.csv" )

  
yName="distance" 
x1Name="epoch" 
x2Name="fish"

#-------------------------------------------------------------------------------

source("run_model_exp1_dir1.R")

expDataFrame = subset(fullDataFrame, experiment == 'exp1_dir1', select = c("distance","fish","epoch"))

mcmcCoda = runMCMC( datFrm=expDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=10000 , thinSteps=1 , nChains = 4)


mcmcMat = as.matrix(mcmcCoda,chains=TRUE)


write.csv(mcmcMat, "post_epoch_subject_dir1.csv")

# diagMCMC( mcmcCoda , parName="b0" )
# diagMCMC( mcmcCoda , parName="bE[1]" )
# diagMCMC( mcmcCoda , parName="bS[1]" )

#-------------------------------------------------------------------------------

source("run_model_exp1_dir2.R")

expDataFrame = subset(fullDataFrame, experiment == 'exp1_dir2', select = c("distance","fish","epoch"))

mcmcCoda = runMCMC( datFrm=expDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=10000 , thinSteps=1 , nChains = 4)

mcmcMat = as.matrix(mcmcCoda,chains=TRUE)

write.csv(mcmcMat, "post_epoch_subject_dir2.csv")

# diagMCMC( mcmcCoda , parName="b0" )
# diagMCMC( mcmcCoda , parName="bE[1]" )
# diagMCMC( mcmcCoda , parName="bS[1]" )

#--------------------------------------------------------------------------------

source("run_model_exp2_reverse.R")

expDataFrame = subset(fullDataFrame, experiment == 'exp2_reverse', select = c("distance","fish","epoch"))

mcmcCoda = runMCMC( datFrm=expDataFrame , 
                    yName=yName , x1Name=x1Name , x2Name=x2Name ,
                    numSavedSteps=10000 , thinSteps=1 , nChains = 4)

mcmcMat = as.matrix(mcmcCoda,chains=TRUE)

write.csv(mcmcMat, "post_epoch_subject_reverse.csv")

# diagMCMC( mcmcCoda , parName="b0" )
# diagMCMC( mcmcCoda , parName="bE[1]" )
# diagMCMC( mcmcCoda , parName="bS[1]" )
