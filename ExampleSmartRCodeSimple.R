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
DataWideFormat$WeightsUsed <- DataWideFormat$KnownWeight1 * DataWideFormat$KnownWeight2;  

################################################################
# Translate wide-format dataset into a long-format dataset;
################################################################
nwaves <- 3; 
DataLongFormat <- reshape(DataWideFormat, 
                          varying = c("Y1", "Y2", "Y3"), 
                          v.names = "Y",
                          timevar = "time", 
                          times = c(1, 2, 3), 
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
DataForAnalysis$S1 <- 1*(DataForAnalysis$time>1);  
# that is, if time > 1 then S1 <- 1; otherwise S1 = 0;
DataForAnalysis$S2 <- 1*(DataForAnalysis$time>2);  
# that is, if time > 2 then S2 <- 1; otherwise S2 = 0;
# Note that time = S1 + S2 + 1;


################################################################
# Do a GEE analysis approach with independent working correlation structure:
################################################################
formula1 <- Y ~ S1 + S2 + S1:A1 +  S2:A1 + S2:A2 + S2:A1:A2;
GeeFinal <- geeglm(formula = formula1,  
                         family=binomial,
                         id=id,
                         weights = WeightsUsed,    
                         data=DataForAnalysis,
                         corstr = "independence"); 
beta <- coef(GeeFinal);  # These are the GEE regression coefficient estimates;
CovBeta <- GeeFinal$geese$vbeta; 
    # This is the covariance matrix of the GEE regression coefficient estimates
    # as provided by the package, but it is not accurate unless known weights
    # are used.;
 

################################################################
# Now that we have beta and CovBeta, we can calculate the estimands of interest,
# such as contrasts between areas under the curve.
################################################################
L.plus.plus   <- rbind(c( 1, 0, 0, 0, 0, 0, 0),
                       c( 1, 1, 0, 1, 0, 0, 0),
                       c( 1, 1, 1, 1, 1, 1, 1));
  # This three-by-seven matrix gives the contrast coefficients (multipliers) for
  # the linear combinations of the seven GEE regression coefficients which will
  # give the inverse logit of the expected probability of Y=1 at each of the 
  # three measurement time points, for the A1=+1, A2=+1 embedded adaptive intervention.
L.plus.minus  <- rbind(c( 1, 0, 0, 0, 0, 0, 0),
                       c( 1, 1, 0, 1, 0, 0, 0),
                       c( 1, 1, 1, 1, 1,-1,-1));
  # This matrix is for the A1=+1, A2=-1 embedded adaptive intervention.
L.minus.plus  <- rbind(c( 1, 0, 0, 0, 0, 0, 0),
                       c( 1, 1, 0,-1, 0, 0, 0),
                       c( 1, 1, 1,-1,-1, 1,-1));
  # This matrix is for the A1=-1, A2=+1 embedded adaptive intervention.
L.minus.minus <- rbind(c( 1, 0, 0, 0, 0, 0, 0),
                       c( 1, 1, 0,-1, 0, 0, 0),
                       c( 1, 1, 1,-1,-1,-1, 1));
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
estimated.area.plus.plus <- sum(c(0.5,1.0,0.5)*mu.plus.plus);
  # Estimated area under the curve for the A1=+1, A2=+1 embedded adaptive intervention.
estimated.area.plus.minus <- sum(c(0.5,1.0,0.5)*mu.plus.minus);
  # Estimated area under the curve for the A1=+1, A2=-1 embedded adaptive intervention.
estimated.area.minus.plus <- sum(c(0.5,1.0,0.5)*mu.minus.plus);
  # Estimated area under the curve for the A1=-1, A2=+1 embedded adaptive intervention.
estimated.area.minus.minus <- sum(c(0.5,1.0,0.5)*mu.minus.minus);
  # Estimated area under the curve for the A1=-1, A2=-1 embedded adaptive intervention.
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
derivatives.pp.versus.pm <- c(.5,1,.5,-.5,-1,-.5)*c(mu.plus.plus, mu.plus.minus)*(1-c(mu.plus.plus, mu.plus.minus));
  # These are used in Cramer's delta method for converting the covariance in the 
  # GEE coefficients to the variance of a contrast of interest.;
derivatives.pp.versus.mp <- c(.5,1,.5,-.5,-1,-.5)*c(mu.plus.plus, mu.minus.plus)*(1-c(mu.plus.plus, mu.minus.plus));
derivatives.pp.versus.mm <- c(.5,1,.5,-.5,-1,-.5)*c(mu.plus.plus, mu.minus.minus)*(1-c(mu.plus.plus, mu.minus.minus));
derivatives.pm.versus.mp <- c(.5,1,.5,-.5,-1,-.5)*c(mu.plus.minus,mu.minus.plus)*(1-c(mu.plus.minus,mu.minus.plus));
derivatives.pm.versus.mm <- c(.5,1,.5,-.5,-1,-.5)*c(mu.plus.minus,mu.minus.minus)*(1-c(mu.plus.minus,mu.minus.minus));
derivatives.mp.versus.mm <- c(.5,1,.5,-.5,-1,-.5)*c(mu.minus.plus,mu.minus.minus)*(1-c(mu.minus.plus,mu.minus.minus));
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
 