################################################################
# Load the geepack library, which you will need:
################################################################
library(geepack);


################################################################
# Read the sample dataset:
################################################################
DataWideFormat <- read.table("SimulatedSmartBinaryData.txt", header=TRUE, na=".");


################################################################
# Create the weights based on known assignment probabilities;
################################################################
DataWideFormat$KnownWeight1 <- 2;
DataWideFormat$KnownWeight2 <- 2*(DataWideFormat$R==0)+1*(DataWideFormat$R==1);  
# that is, IF (R = 0) THEN KnownWeight2 = 2; ELSE KnownWeight2 = 1;
DataWideFormat$KnownWeight <- DataWideFormat$KnownWeight1 * DataWideFormat$KnownWeight2; 


################################################################
# If desired, also create weights based on estimated assignment probabilities.
################################################################
# Estimate Stage 1 weight:
DataWideFormat$A1DummyCoded <- 0*(DataWideFormat$A1==-1)+1*(DataWideFormat$A1==+1);
logisticModel1 <- glm(formula=A1DummyCoded ~  Male + BaselineSeverity,    
                           # You could include baseline covariates in this formula if 
                           # they may be relevant.
					family=binomial,
					data=DataWideFormat);
DataWideFormat$p1 <- fitted(logisticModel1);
DataWideFormat$EstimatedWeight1 <- DataWideFormat$A1DummyCoded/DataWideFormat$p1 +
(1-DataWideFormat$A1DummyCoded)/(1-DataWideFormat$p1);
# Estimate Stage 2 weight:
DataWideFormat$A2DummyCoded <- 0*(DataWideFormat$A2==-1)+1*(DataWideFormat$A2==+1);
DataWideFormat$A2DummyCoded[which(DataWideFormat$R==1)] <- NA;
# responders were not re-randomized, so A2 is neither -1 nor +1;
logisticModel2 <- glm(formula=A2DummyCoded ~  Y1 + Male + BaselineSeverity,  
					           # You could include baseline covariates or other aspects of treatment history
					           # in this formula if they may be relevant.   However, you cannot include 
					           # response status R as a predictor of A2, because A2 only has variance if R=0.
					family=binomial,
					na.action=na.exclude,
					data=DataWideFormat); 
DataWideFormat$p2 <- fitted(logisticModel2);
DataWideFormat$EstimatedWeight2 <- 1;     
who.responded <- which(DataWideFormat$R==0);
DataWideFormat$EstimatedWeight2[who.responded] <- 
DataWideFormat$A2DummyCoded[who.responded]/DataWideFormat$p2[who.responded] +
      (1-DataWideFormat$A2DummyCoded[who.responded])/(1-DataWideFormat$p2[who.responded]);
# Responders have a stage 2 weight of 1;  nonresponders have a stage 2 weighted which is their
# estimated probability of getting the treatment they did in fact get.
DataWideFormat$EstimatedWeight <- DataWideFormat$EstimatedWeight1 * DataWideFormat$EstimatedWeight2;
# Note:  Because the estimated weights are not necessarily integers, you may get a warning message
# from R when using them in a logistic regression :  "non-integer #successes in a binomial glm!"
# You can ignore this.

################################################################
# Choose which weights you want to use.
################################################################
# For this part, either say 
#    DataWideFormat$WeightsUsed <- DataWideFormat$KnownWeight;
# or 
#    DataWideFormat$WeightsUsed <- DataWideFormat$EstimatedWeight;
DataWideFormat$WeightsUsed <- DataWideFormat$EstimatedWeight;


################################################################
# Translate wide-format dataset into a long-format dataset;
################################################################
nwaves <- 6; 
DataLongFormat <- reshape(DataWideFormat, 
                          varying = c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6"), 
                          v.names = "Y",
                          timevar = "time", 
                          times = 1:6, 
                          new.row.names = 1:(nwaves*nrow(DataWideFormat)),
                          direction = "long");


################################################################
# Create the "replications."
################################################################
# That is, replicate people who got A1=+1 and responded with R=1.  They were not
# re-randomized, so we have to replicate their data to inform both of
# the dynamic treatment regimens which they could have received had they
# been re-randomized.  It is somewhat as if we are creating two clones of
# each of them and counting each in a different treatment, because in
# reality it is indistinguishable which one they received.
RowsToReplicate <- DataLongFormat[which(DataLongFormat$R==1),];
RowsNotToReplicate <- DataLongFormat[which(DataLongFormat$R==0),];
RowsNotToReplicate$wave <- RowsNotToReplicate$time;  
# "wave" will be a recoded version of "time" which is different for
# some of the replicated pseudodata.
PlusOnePseudodata <- RowsToReplicate;
PlusOnePseudodata$A2 <- 1;
PlusOnePseudodata$wave <- PlusOnePseudodata$time;
MinusOnePseudodata <- RowsToReplicate;
MinusOnePseudodata$A2 <- -1;
MinusOnePseudodata$wave <- MinusOnePseudodata$time + nwaves;  
# We keep the same subject ID to show that we don't really have all those
# new participants.  So we have to distinguish the new observations somehow,
# and so we treat them as new waves of data on the same person.  Although
# it seems very ad-hoc, this method has been shown to be valid.
# Create the final analysis dataset including replicates.
DataForAnalysis <- rbind(PlusOnePseudodata, MinusOnePseudodata, RowsNotToReplicate);
DataForAnalysis <- DataForAnalysis[order(DataForAnalysis$id,DataForAnalysis$wave),] 
# This sorts by the variables id and wave
DataForAnalysis$S1 <- 0.5 + 1*(DataForAnalysis$time>1);  # time since first randomization;
# that is, if time > 1 then S1 <- 1.5; otherwise S1 = 0.5;
DataForAnalysis$S2 <- pmax(0,DataForAnalysis$time-2);   # time since second randomization;
# that is, if time > 2 then S2 <- time-2; otherwise S2 = 0; 


################################################################
# Do a GEE analysis approach with independent working correlation structure:
################################################################
GeeModelFormula <-  Y ~ Male + BaselineSeverity +
  S1 + S2 + S1:A1 + S2:A1 + S2:A2 + 
  S2:A1:A2 ;
GeeIndependent <- geeglm(formula = GeeModelFormula,  
                         family=binomial,
                         id=id,
                         weights = WeightsUsed,    
                         data=DataForAnalysis,
                         corstr = "independence");


################################################################
# You can skip this part if you just want to use working independence.
################################################################
# You may be able to improve upon the working-independence results by
# doing a GEE analysis approach with block-exchangeable or block-AR-1
# correlation structure.
# The approach we use here involves an initial fit of an independence 
# model to get a correlation parameter estimate, which is then used
# in a second model fit.  
# The code for this approach assumes that there are no missing observations, 
# because it involves lining up observations from different waves to estimate 
# covariance.  The independence approach above does not need this assumption.  
# A more complicated version of this code could remove this assumption, 
# but we have tried to make this tutorial code simple. 
# The first part of the approach is the same whether you want to use AR-1 or exchangeable.
# We create a table of estimated variances indexed by A1, A2, and time.;
ResidualVarianceByRegimenAndTime <- expand.grid(A1=unique(DataForAnalysis$A1),
                                                A2=unique(DataForAnalysis$A2),
                                                time=unique(DataForAnalysis$time),
                                                VarianceEstimate=NA);
for (i in 1:nrow(ResidualVarianceByRegimenAndTime)) {
  TheseRows <- which( (DataForAnalysis$A1==ResidualVarianceByRegimenAndTime$A1[i])&
                        (DataForAnalysis$A2==ResidualVarianceByRegimenAndTime$A2[i])&
                        (DataForAnalysis$time==ResidualVarianceByRegimenAndTime$time[i]) );   
  ResidualVarianceByRegimenAndTime$VarianceEstimate[i] <- 
    weighted.mean(x=(GeeIndependent$resid[TheseRows]^2),w=DataForAnalysis$WeightsUsed[TheseRows]);
}  
# Now we create a smaller table of estimated variances indexed only by A1 and A2.
ResidualCovarianceByRegimen <- expand.grid(A1=unique(DataForAnalysis$A1),
                                           A2=unique(DataForAnalysis$A2), 
                                           VarianceEstimate=NA,
                                           CrossCorrelationEstimate=NA);
for (i in 1:nrow(ResidualCovarianceByRegimen)) {
  TheseRows <- which( (ResidualVarianceByRegimenAndTime$A1==ResidualCovarianceByRegimen$A1[i])&
                        (ResidualVarianceByRegimenAndTime$A2==ResidualCovarianceByRegimen$A2[i])  );   
  ResidualCovarianceByRegimen$VarianceEstimate[i] <- 
    mean(ResidualVarianceByRegimenAndTime$VarianceEstimate[TheseRows]);
}   
VarianceEstimatePooled <- mean(ResidualCovarianceByRegimen$VarianceEstimate);


################################################################
# Now you could continue with this code if you want to use exchangeable working correlation...
################################################################
# Estimate off-diagonal within-person correlation:
for (i in 1:nrow(ResidualCovarianceByRegimen)) {  # regimen index;
  MeanResidualCrossProductsPerSubjectForThisRegimen <- NULL;
  WeightsForThisRegimen <- NULL;
  ThisRegimen <- which( (DataForAnalysis$A1==ResidualCovarianceByRegimen$A1[i])&
                          (DataForAnalysis$A2==ResidualCovarianceByRegimen$A2[i]));
  # Calculate the exchangeable correlation parameter (if you want to use exchangeable correlation): 
  for (j in unique(DataForAnalysis$id)) { # subject index; 
    ThisSubjectAndRegimen <- which( (DataForAnalysis$A1==ResidualCovarianceByRegimen$A1[i])&
                                      (DataForAnalysis$A2==ResidualCovarianceByRegimen$A2[i]) &
                                      (DataForAnalysis$id==j));
    if (length(ThisSubjectAndRegimen)>0) {
      tempCrossProducts <- crossprod(t(GeeIndependent$resid[ThisSubjectAndRegimen]))[
        upper.tri(crossprod(t(GeeIndependent$resid[ThisSubjectAndRegimen])))]; 
      MeanResidualCrossProductsPerSubjectForThisRegimen <- c(
                               MeanResidualCrossProductsPerSubjectForThisRegimen,
                               mean(tempCrossProducts)); 
      WeightsForThisRegimen <- c(WeightsForThisRegimen, DataForAnalysis[ThisSubjectAndRegimen[1],]$WeightsUsed);
      # assumes all observations per subject have the same weight;
    }       
  } 
  ResidualCovarianceByRegimen$CrossCorrelationEstimate[i] <- weighted.mean(
                             x=MeanResidualCrossProductsPerSubjectForThisRegimen,
                             w=WeightsForThisRegimen)/ 
    ResidualCovarianceByRegimen$VarianceEstimate[i];
} 
CrossCorrelationEstimatePooled <- mean(ResidualCovarianceByRegimen$CrossCorrelationEstimate);
rho <- CrossCorrelationEstimatePooled;
ExchangeableWorkingCorrelation <- diag(as.vector(rep(1-rho,nwaves)))+matrix(rho,nwaves,nwaves);
BlockStructureWorkingCorrelation <- rbind(cbind(ExchangeableWorkingCorrelation,
                                              0*ExchangeableWorkingCorrelation),
                                        cbind(0*ExchangeableWorkingCorrelation,
                                              ExchangeableWorkingCorrelation)); 
WorkingCorrelationInGeepackFormat <- fixed2Zcor(cor.fixed=BlockStructureWorkingCorrelation, 
                                                id=DataForAnalysis$id,  
                                                waves=DataForAnalysis$wave); 
GeeExchangeable <- geeglm(formula = GeeModelFormula,  
                          id=id,  
                          family=binomial,
                          weights = WeightsUsed,
                          data=DataForAnalysis,
                          corstr = "fixed",
                          zcor=WorkingCorrelationInGeepackFormat); 


################################################################
# Or you could continue with this code if you want to use AR-1 working correlation.
# It is similar to the working-exchangeable code above.
################################################################
# Calculate the AR-1 parameter (if you are using AR-1) 
for (j in unique(DataForAnalysis$id)) { # subject index; 
  ThisSubjectAndRegimen <- which( (DataForAnalysis$A1==ResidualCovarianceByRegimen$A1[i])&
                                    (DataForAnalysis$A2==ResidualCovarianceByRegimen$A2[i]) &
                                    (DataForAnalysis$id==j));
  if (length(ThisSubjectAndRegimen)>0) {
    tempVector <- GeeIndependent$resid[ThisSubjectAndRegimen];
    tempCrossProducts <- tempVector[1:(length(tempVector)-1)] * tempVector[2:length(tempVector)]; 
    MeanResidualCrossProductsPerSubjectForThisRegimen <- c(MeanResidualCrossProductsPerSubjectForThisRegimen,
                                                           mean(tempCrossProducts)); 
    WeightsForThisRegimen <- c(WeightsForThisRegimen, DataForAnalysis[ThisSubjectAndRegimen[1],]$WeightsUsed);
    # assumes all observations per subject have the same weight;
  }       
} 
ResidualCovarianceByRegimen$CrossCorrelationEstimate[i] <- weighted.mean(x=MeanResidualCrossProductsPerSubjectForThisRegimen,
                                                                         w=WeightsForThisRegimen)/ 
  ResidualCovarianceByRegimen$VarianceEstimate[i]; 
CrossCorrelationEstimatePooled <- mean(ResidualCovarianceByRegimen$CrossCorrelationEstimate); 
rho <- CrossCorrelationEstimatePooled;  
AR1WorkingCorrelation <- matrix(0,nwaves,nwaves);
for (thisRow in 1:nrow(AR1WorkingCorrelation)) {
  for (thisColumn in 1:ncol(AR1WorkingCorrelation)) {
    AR1WorkingCorrelation[thisRow,thisColumn] <- rho^abs(thisRow-thisColumn);
  }
} 
BlockStructureWorkingCorrelation <- rbind(cbind(AR1WorkingCorrelation,0*AR1WorkingCorrelation),
                               cbind(0*AR1WorkingCorrelation,AR1WorkingCorrelation)); 
WorkingCorrelationInGeepackFormat <- fixed2Zcor(cor.fixed=BlockStructureWorkingCorrelation, 
                                                id=DataForAnalysis$id,  
                                                waves=DataForAnalysis$wave); 
GeeAR1 <- geeglm(formula = GeeModelFormula,  
                 id=id,  
                 family=binomial,
                 weights = WeightsUsed,
                 data=DataForAnalysis,
                 corstr = "fixed",
                 zcor=WorkingCorrelationInGeepackFormat); 


################################################################
# Regardless of the working structure you choose to use, you will
# have to do some more steps to get the estimands of interest and their
# standard error.  Those steps require the GEE result, which 
# we will call GeeFinal here.  You could
# choose 
#    GeeFinal <- GeeIndependent;
#    GeeFinal <- GeeExchangeable;
#    or
#    GeeFinal <- GeeAR1;
################################################################
GeeFinal <- GeeAR1;
beta <- coef(GeeFinal);  # These are the GEE regression coefficient estimates;
CovBeta <- GeeFinal$geese$vbeta; 
    # This is the covariance matrix of the GEE regression coefficient estimates
    # as provided by the package, but it is not accurate unless known weights
    # are used.;


################################################################
# If you are using estimated weights, recalculate CovBeta taking into account the  
# use of the estimated weights.  Skip this section if you are not using 
# estimated weights.;
################################################################
nsub <- nrow(DataWideFormat);
PredictorsLogistic1 =  cbind(DataWideFormat$Male, DataWideFormat$BaselineSeverity); 
PredictorsLogistic2 = cbind( DataWideFormat$Y1,
                             DataWideFormat$Male, DataWideFormat$BaselineSeverity); 
scoreAlpha1 <- cbind(1,PredictorsLogistic1) * (DataWideFormat$A1DummyCoded - DataWideFormat$p1);
residuals2 <- DataWideFormat$A2DummyCoded - DataWideFormat$p2;
residuals2[which(is.na(residuals2))] <- 0;  # This excludes the responders from the model for stage 2;
scoreAlpha2 <- cbind(1,PredictorsLogistic2) * residuals2; 
# Step 3: Finish calculating the standard errors.
# This is done with the replicated data. 
predictors <- cbind(DataForAnalysis$Male,
                    DataForAnalysis$BaselineSeverity,
                    DataForAnalysis$S1,
                    DataForAnalysis$S2,
                    DataForAnalysis$S1*DataForAnalysis$A1,
                    DataForAnalysis$S2*DataForAnalysis$A1,
                    DataForAnalysis$S2*DataForAnalysis$A2,
                    DataForAnalysis$S2*DataForAnalysis$A1*DataForAnalysis$A2);
predictors <- cbind(1,predictors); # add intercept column;
# Construct "matrixI," the empirical covariance matrix of the 
# GEE coefficients, and "matrixJ", their model-based covariance matrix. 
# "I" and "J" are the notation used in Xi Lu's code (I stands for "information" as in
# Fisher information matrix.  They don't represent the usual I and J matrices which 
#	are the identity (diagonal) and all-ones matrix respectively. */
resid <- GeeFinal$residuals; 
fitted <- GeeFinal$fitted.values; 
stopifnot(length(resid)==nrow(predictors)); 
             # Check to make sure there is a residual for each observation 
             # and that nothing weird happened because of NA's.;
matrixJ <- matrix(0,ncol(predictors),ncol(predictors));
U <- matrix(0,nsub,ncol(predictors));
# The rows of U are the score function for each subject (derivative of that subject's
# log-pseudo-likelihood contribution with respect to beta), with each subject in the original 
# data set having one row. */   
indicesNextSubject <- 1:nwaves;  
# In the upcoming look, you will need the estimated nwaves*nwaves correlation matrix of the data.
# If you chose to use working independence, use
#    WorkCorr <- diag(nwaves);
# If you chose to use working exchangeable, use 
#    WorkCorr <- ExchangeableWorkingCorrelation;
# If you chose to use working autoregressive AR-1, use 
#    WorkCorr <- AR1WorkingCorrelation;
WorkCorr <- AR1WorkingCorrelation;
for (i in 1:nsub) {
  indicesThisSubject <- indicesNextSubject;
  weightThisSubject <- DataForAnalysis$KnownWeight[indicesThisSubject[1]];
  # The weights are the same for all observations within the 
  # subject, so we just read the first one. 
  residualsThisSubject <- resid[indicesThisSubject];
  predictorsThisSubject <- predictors[indicesThisSubject,] ;  
  fittedValuesThisSubject <- fitted[indicesThisSubject];
  Mi <- fittedValuesThisSubject*(1-fittedValuesThisSubject); 
  errorThisSubject <- weightThisSubject * t(predictorsThisSubject) %*% 
    diag(as.vector(sqrt(Mi)))%*% 
    solve(WorkCorr) %*% 
    diag(as.vector(1/sqrt(Mi)))%*%
    residualsThisSubject;
  # errorThisSubject is a k by 1 vector for this subject's score function, 
  # here k is ncol(predictors) including intercept. ;
  if (indicesThisSubject[nwaves] == nrow(DataForAnalysis)) {
    # This is the case where the current observations are for the last individual,
    # so no more observations are forthcoming. 
    U[i,] <- t(errorThisSubject);
    thisMatrixJ <- weightThisSubject * t(predictorsThisSubject)%*%
      diag(as.vector(sqrt(Mi)))%*% 
      solve(WorkCorr) %*% 
      diag(as.vector(1/sqrt(Mi)))%*%
      predictorsThisSubject;
    # contribution to J from this subject;
    matrixJ <- matrixJ + thisMatrixJ;
  } else {
    if (DataForAnalysis$wave[[indicesThisSubject[nwaves]+1]] == 1) {
      # This is the case where the next observations are from an actual new individual. */
      U[i,] <- t(errorThisSubject);
      thisMatrixJ = weightThisSubject * t(predictorsThisSubject)%*%
        diag(as.vector(sqrt(Mi)))%*% 
        solve(WorkCorr) %*% 
        diag(as.vector(1/sqrt(Mi)))%*%
        predictorsThisSubject;
      # contribution to J from this subject;
      matrixJ <- matrixJ + thisMatrixJ;
      indicesNextSubject <- (indicesThisSubject[nwaves]+1):(indicesThisSubject[nwaves]+nwaves); 
    } 
    if (DataForAnalysis$wave[[indicesThisSubject[nwaves]+1]] == nwaves+1) {
      # This is the case where the next observations are a replicate of this individual.;
      indicesReplicate <- (indicesThisSubject[nwaves]+1):(indicesThisSubject[nwaves]+nwaves);
      residReplicate <- resid[indicesReplicate];
      predictorsReplicate <- predictors[indicesReplicate,];
      errorReplicate <- weightThisSubject * t(predictorsReplicate)%*%
        diag(as.vector(sqrt(Mi)))%*% 
        solve(WorkCorr) %*% 
        diag(as.vector(1/sqrt(Mi)))%*%
        residReplicate;
      U[i,] <- t(errorThisSubject + errorReplicate);
      thisMatrixJ <- weightThisSubject * t(predictorsThisSubject)%*%
        diag(as.vector(sqrt(Mi)))%*% 
        solve(WorkCorr) %*% 
        diag(as.vector(1/sqrt(Mi)))%*%
        predictorsThisSubject + 
        weightThisSubject * t(predictorsReplicate)%*%
        diag(as.vector(sqrt(Mi)))%*% 
        solve(WorkCorr) %*% 
        diag(as.vector(1/sqrt(Mi)))%*%
        predictorsReplicate;
      matrixJ <- matrixJ + thisMatrixJ;
      indicesNextSubject <- (indicesReplicate[nwaves]+1):(indicesReplicate[nwaves]+nwaves);
    }
    if ((DataForAnalysis$wave[indicesThisSubject[nwaves]+1] != 1) & 
        (DataForAnalysis$wave[indicesThisSubject[nwaves]+1] != nwaves+1) ) {
      print("Unexpected error in code or dataset;");
      print(indicesThisSubject); 
    }
  }
}
#  Construct the final covariance matrix estimate using an adjusted sandwich formula, 
# and extract the standard errors.;
matrixIKnownWeights <- (1/nsub)*t(U)%*%U;
scores <- cbind(scoreAlpha1, scoreAlpha2);
hat <- scores%*%solve(t(scores)%*%scores)%*%t(scores)
# The "hat" or projection matrix of the regression of the weights on the predictors of the weights.
Uproj <- U - hat%*%U;
matrixI <- (1/nsub)*t(Uproj)%*%Uproj;
matrixIAdjusted <- (nsub/(nsub-ncol(predictors)))*matrixI; 
    # Adjust variance estimate for presence of covariates;
matrixJ <- matrixJ/nsub;
invJ <- solve(matrixJ); 
CovBeta <- (1/nsub)* invJ %*% matrixIAdjusted %*% invJ;

################################################################
# Now that we have beta and CovBeta, we can calculate the estimands of interest,
# such as contrasts between areas under the curve.
################################################################
L.plus.plus   <- matrix(c(
  1, 1, 1, 0.5, 0, +0.5,  0,  0,  0,
  1, 1, 1, 1.5, 0, +1.5,  0,  0,  0,
  1, 1, 1, 1.5, 1, +1.5, +1, +1, +1,
  1, 1, 1, 1.5, 2, +1.5, +2, +2, +2,
  1, 1, 1, 1.5, 3, +1.5, +3, +3, +3,
  1, 1, 1, 1.5, 4, +1.5, +4, +4, +4), nrow=6, byrow=TRUE);
  # This six-by-nine matrix gives the contrast coefficients (multipliers) for
  # the linear combinations of the nine GEE regression coefficients which will
  # give the inverse logit of the expected probability of Y=1 at each of the 
  # six measurement time points, for the A1=+1, A2=+1 embedded adaptive intervention.
L.plus.minus  <- matrix(c(
  1, 1, 1, 0.5, 0, +0.5,  0,  0,  0,
  1, 1, 1, 1.5, 0, +1.5,  0,  0,  0,
  1, 1, 1, 1.5, 1, +1.5, +1, -1, -1,
  1, 1, 1, 1.5, 2, +1.5, +2, -2, -2,
  1, 1, 1, 1.5, 3, +1.5, +3, -3, -3,
  1, 1, 1, 1.5, 4, +1.5, +4, -4, -4),nrow=6,byrow=TRUE);
  # This matrix is for the A1=+1, A2=-1 embedded adaptive intervention.
L.minus.plus  <- matrix(c(
  1, 1, 1, 0.5, 0, -0.5,  0,  0,  0,
  1, 1, 1, 1.5, 0, -1.5,  0,  0,  0,
  1, 1, 1, 1.5, 1, -1.5, -1, +1, -1,
  1, 1, 1, 1.5, 2, -1.5, -2, +2, -2,
  1, 1, 1, 1.5, 3, -1.5, -3, +3, -3,
  1, 1, 1, 1.5, 4, -1.5, -4, +4, -4),nrow=6,byrow=TRUE);
  # This matrix is for the A1=-1, A2=+1 embedded adaptive intervention.
L.minus.minus <- matrix(c(
  1, 1, 1, 0.5, 0, -0.5,  0,  0,  0,
  1, 1, 1, 1.5, 0, -1.5,  0,  0,  0,
  1, 1, 1, 1.5, 1, -1.5, -1, -1, +1,
  1, 1, 1, 1.5, 2, -1.5, -2, -2, +2,
  1, 1, 1, 1.5, 3, -1.5, -3, -3, +3,
  1, 1, 1, 1.5, 4, -1.5, -4, -4, +4), nrow=6, byrow=TRUE);
  # This matrix is for the A1=-1, A2=-1 embedded adaptive intervention.
eta.plus.plus <- as.vector(L.plus.plus%*%beta);  
  # vector of fitted inverse-logit probabilities of Y=1 at each time point for (+1,+1);
eta.plus.minus <- as.vector(L.plus.minus%*%beta);
  # vector of fitted inverse-logit probabilities of Y=1 at each time point for (+1,-1);
eta.minus.plus <- as.vector(L.minus.plus%*%beta);
  # vector of fitted inverse-logit probabilities of Y=1 at each time point for (-1,+1);
eta.minus.minus <- as.vector(L.minus.minus%*%beta);
  # vector of fitted inverse-logit probabilities of Y=1 at each time point for (-1,-1);
mu.plus.plus <- exp(eta.plus.plus)/(1+exp(eta.plus.plus));
  # vector of fitted probabilities of Y=1 at each time point for (+1,+1);
mu.plus.minus <- exp(eta.plus.minus)/(1+exp(eta.plus.minus));
  # vector of fitted probabilities of Y=1 at each time point for (+1,-1);
mu.minus.plus <- exp(eta.minus.plus)/(1+exp(eta.minus.plus));
  # vector of fitted probabilities of Y=1 at each time point for (-1,+1);
mu.minus.minus <- exp(eta.minus.minus)/(1+exp(eta.minus.minus));
  # vector of fitted probabilities of Y=1 at each time point for (-1,-1);
estimated.area.plus.plus <- sum((1/5)*c(0.5,1,1,1,1,0.5)*mu.plus.plus);
# Estimated area under the curve for the A1=+1, A2=+1 embedded adaptive intervention,
# divided by the length of the time interval.
estimated.area.plus.minus <- sum((1/5)*c(0.5,1,1,1,1,0.5)*mu.plus.minus);
# Estimated area under the curve for the A1=+1, A2=-1 embedded adaptive intervention,
# divided by the length of the time interval.
estimated.area.minus.plus <- sum((1/5)*c(0.5,1,1,1,1,0.5)*mu.minus.plus);
# Estimated area under the curve for the A1=-1, A2=+1 embedded adaptive intervention,
# divided by the length of the time interval.
estimated.area.minus.minus <- sum((1/5)*c(0.5,1,1,1,1,0.5)*mu.minus.minus);
# Estimated area under the curve for the A1=-1, A2=-1 embedded adaptive intervention,
# divided by the length of the time interval.
estimate.pp.versus.pm <- estimated.area.plus.plus-estimated.area.plus.minus;
  # The estimated difference equals the difference in the estimates;
estimate.pp.versus.mp <- estimated.area.plus.plus-estimated.area.minus.plus;
estimate.pp.versus.mm <- estimated.area.plus.plus-estimated.area.minus.minus;
estimate.pm.versus.mp <- estimated.area.plus.minus-estimated.area.minus.plus;
estimate.pm.versus.mm <- estimated.area.plus.minus-estimated.area.minus.minus;
estimate.mp.versus.mm <- estimated.area.minus.plus-estimated.area.minus.minus;
estimated.values.for.contrasts <- c(estimate.pp.versus.pm=estimate.pp.versus.pm,
                                    estimate.pp.versus.mp=estimate.pp.versus.mp,
                                    estimate.pp.versus.mm=estimate.pp.versus.mm,
                                    estimate.pm.versus.mp=estimate.pm.versus.mp,
                                    estimate.pm.versus.mm=estimate.pm.versus.mm,
                                    estimate.mp.versus.mm=estimate.mp.versus.mm); 
  # These are the choose(4,2) = 6 contrasts of interest among the AUC's for the
  # four embedded adaptive interventions.;
names(estimated.values.for.contrasts) <- c("AUC, ++ versus +-",
                                           "AUC, ++ versus -+",
                                           "AUC, ++ versus --",
                                           "AUC, +- versus -+",
                                           "AUC, +- versus --",
                                           "AUC, -+ versus --");
derivatives.pp.versus.pm <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.plus.plus, mu.plus.minus)*(1-c(mu.plus.plus, mu.plus.minus));
  # These are used in Cramer's delta method for converting the covariance in the 
  # GEE coefficients to the variance of a contrast of interest.;
derivatives.pp.versus.mp <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.plus.plus, mu.minus.plus)*(1-c(mu.plus.plus, mu.minus.plus));
derivatives.pp.versus.mm <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.plus.plus, mu.minus.minus)*(1-c(mu.plus.plus, mu.minus.minus));
derivatives.pm.versus.mp <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.plus.minus,mu.minus.plus)*(1-c(mu.plus.minus,mu.minus.plus));
derivatives.pm.versus.mm <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.plus.minus,mu.minus.minus)*(1-c(mu.plus.minus,mu.minus.minus));
derivatives.mp.versus.mm <- (1/5)*c(.5,1,1,1,1,.5,-.5,-1,-1,-1,-1,-.5)*c(mu.minus.plus,mu.minus.minus)*(1-c(mu.minus.plus,mu.minus.minus));
multipliers.pp.versus.pm <- rbind(L.plus.plus, L.plus.minus);
multipliers.pp.versus.mp <- rbind(L.plus.plus, L.minus.plus);
multipliers.pp.versus.mm <- rbind(L.plus.plus, L.minus.minus);
multipliers.pm.versus.mp <- rbind(L.plus.minus,L.minus.plus);
multipliers.pm.versus.mm <- rbind(L.plus.minus,L.minus.minus);
multipliers.mp.versus.mm <- rbind(L.minus.plus,L.minus.minus);
stderr.pp.versus.pm <- sqrt(derivatives.pp.versus.pm%*%
                              multipliers.pp.versus.pm%*%
                              CovBeta%*%
                              t(multipliers.pp.versus.pm)%*%
                              derivatives.pp.versus.pm);
   # This is an example of Cramer's delta method;
stderr.pp.versus.mp <- sqrt(derivatives.pp.versus.mp%*%
                              multipliers.pp.versus.mp%*%
                              CovBeta%*%
                              t(multipliers.pp.versus.mp)%*%
                              derivatives.pp.versus.mp);
stderr.pp.versus.mm <- sqrt(derivatives.pp.versus.mm%*%
                              multipliers.pp.versus.mm%*%
                              CovBeta%*%
                              t(multipliers.pp.versus.mm)%*%
                              derivatives.pp.versus.mm);
stderr.pm.versus.mp <- sqrt(derivatives.pm.versus.mp%*%
                              multipliers.pm.versus.mp%*%
                              CovBeta%*%
                              t(multipliers.pm.versus.mp)%*%
                              derivatives.pm.versus.mp);
stderr.pm.versus.mm <- sqrt(derivatives.pm.versus.mm%*%
                              multipliers.pm.versus.mm%*%
                              CovBeta%*%
                              t(multipliers.pm.versus.mm)%*%
                              derivatives.pm.versus.mm);
stderr.mp.versus.mm <- sqrt(derivatives.mp.versus.mm%*%
                              multipliers.mp.versus.mm%*%
                              CovBeta%*%
                              t(multipliers.mp.versus.mm)%*%
                              derivatives.mp.versus.mm);
standard.errors.for.contrasts <- c(stderr.pp.versus.pm=stderr.pp.versus.pm,
                                   stderr.pp.versus.mp=stderr.pp.versus.mp,
                                   stderr.pp.versus.mm=stderr.pp.versus.mm,
                                   stderr.pm.versus.mp=stderr.pm.versus.mp,
                                   stderr.pm.versus.mm=stderr.pm.versus.mm,
                                   stderr.mp.versus.mm=stderr.mp.versus.mm); 


################################################################
# Present a table of results:
################################################################
z.ratios.for.contrasts <- estimated.values.for.contrasts/standard.errors.for.contrasts;
p.values.for.contrasts <- 2*(1-pnorm(abs(z.ratios.for.contrasts)));
print(round(cbind(estimated.values.for.contrasts,
            standard.errors.for.contrasts,
            z.ratios.for.contrasts,
            p.values.for.contrasts),4));