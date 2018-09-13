/***************************************************************
* Prepare the Dataset
***************************************************************/
* Read the sample dataset;
PROC IMPORT OUT=DataWideFormat 
            DATAFILE= "C:\Users\jjd264\Documents\LongitudinalSmart\BinaryPaper\Demonstration\SimulatedSmartBinaryData.txt" 
            DBMS=TAB REPLACE; 
     GETNAMES=YES;
     DATAROW=2; 
RUN;
* Create known weights;
DATA DataWideFormat; SET DataWideFormat;
	KnownWeight1 = 2;
	IF R = 0 THEN KnownWeight2 = 2; ELSE KnownWeight2 = 1;
	KnownWeight = KnownWeight1 * KnownWeight2; 
RUN;
* Translate wide-format dataset into a long-format dataset;
DATA DataLongFormat;
	SET DataWideFormat;
	Y = Y1;
	time = 1;
	OUTPUT;
	Y = Y2;
	time = 2;
	OUTPUT;
	Y = Y3;
	time = 3;
	OUTPUT;
	Y = Y4;
	time = 4;
	OUTPUT;
	Y = Y5;
	time = 5;
	OUTPUT;
	Y = Y6;
	time = 6;
	OUTPUT;
	KEEP ID KnownWeight Male BaselineSeverity A1 R A2 Y time; 
RUN;
/*****************************************************************************
* Create "replications.";
* That is, replicate people who got A1=+1 and responded with R=1.  They were not;
* re-randomized, so we have to replicate their data to inform both of;
* the dynamic treatment regimens which they could have received had they;
* been re-randomized.  It is somewhat as if we are creating two clones of;
* each of them and counting each in a different treatment, because in;
* reality it is indistinguishable which one they received.;
******************************************************************************/
DATA RowsToReplicate;
	SET DataLongFormat;
	WHERE R=1;
RUN;
DATA RowsNotToReplicate;
	SET DataLongFormat;
	WHERE R=0;
	wave = time;
RUN;
DATA PlusOnePseudodata;
	SET RowsToReplicate;
	A2 = 1;
	wave = time;
RUN;
DATA MinusOnePseudodata;
	SET RowsToReplicate;
	A2 = -1;
	wave = time + 6;
RUN;
DATA DataForAnalysis;
	SET PlusOnePseudodata MinusOnePseudodata RowsNotToReplicate;
RUN;
PROC SORT DATA=DataForAnalysis;
	BY ID Wave;
RUN;
/* Recode time to take changepoint into account: */
DATA DataForAnalysis; SET DataForAnalysis;
	IF time > 1 THEN S1 = 1.5; ELSE S1 = 0.5; /* time since first randomization */
	IF time > 2 THEN S2 = time-2; ELSE S2 = 0;  /* time since second randomization */
RUN;

/***************************************************************
* Do a GEE analysis approach with independent working correlation structure:
***************************************************************/
PROC GENMOD DATA=DataForAnalysis DESCENDING;
	CLASS ID;
	ODS OUTPUT GEEEmpPEst=Beta GEERCov=CovBeta;
	MODEL  Y = Male BaselineSeverity S1 S2 S1*A1 S2*A1 S2*A2 S2*A1*A2 / DIST=BINOMIAL LINK=LOGIT;
	REPEATED SUBJECT=ID / ECOVB;
	WEIGHT KnownWeight;
RUN;

/***************************************************************
* Calculate the estimands of interest and their standard errors. 
***************************************************************/
PROC IML;
	USE Beta;
		READ ALL VAR {Estimate} INTO Beta;
	CLOSE Beta;
	USE CovBeta;
		READ ALL INTO CovBeta;
	CLOSE CovBeta;
	L_Plus_Plus = {   1  1  1  0.5  0  0.5   0   0   0, 
					  1  1  1  1.5  0  1.5   0   0   0, 
					  1  1  1  1.5  1  1.5   1   1   1, 
					  1  1  1  1.5  2  1.5   2   2   2, 
					  1  1  1  1.5  3  1.5   3   3   3, 
					  1  1  1  1.5  4  1.5   4   4   4 };
       * This six-by-nine matrix gives the contrast coefficients (multipliers) for;
	   * the linear combinations of the nine GEE regression coefficients which will;
	   * give the inverse logit of the expected probability of Y=1 at each of the ;
	   * six measurement time points, for the A1=+1, A2=+1 embedded adaptive intervention.;
    L_Plus_Minus = {  1  1  1  0.5  0  0.5   0   0   0, 
					  1  1  1  1.5  0  1.5   0   0   0, 
					  1  1  1  1.5  1  1.5   1  -1  -1, 
					  1  1  1  1.5  2  1.5   2  -2  -2, 
					  1  1  1  1.5  3  1.5   3  -3  -3, 
					  1  1  1  1.5  4  1.5   4  -4  -4 };
       * L_Plus_Minus is for the A1=+1, A2=-1 embedded adaptive intervention.;
	L_Minus_Plus = {  1  1  1  0.5  0  -0.5   0   0   0, 
					  1  1  1  1.5  0  -1.5   0   0   0, 
					  1  1  1  1.5  1  -1.5  -1   1  -1,
					  1  1  1  1.5  2  -1.5  -2   2  -2,
					  1  1  1  1.5  3  -1.5  -3   3  -3, 
					  1  1  1  1.5  4  -1.5  -4   4  -4 };
       * L_Minus_Plus is for the A1=-1, A2=+1 embedded adaptive intervention.;
    L_Minus_Minus = { 1  1  1  0.5  0  -0.5   0   0   0, 
					  1  1  1  1.5  0  -1.5   0   0   0, 
					  1  1  1  1.5  1  -1.5  -1  -1   1, 
					  1  1  1  1.5  2  -1.5  -2  -2   2, 
					  1  1  1  1.5  3  -1.5  -3  -3   3, 
					  1  1  1  1.5  4  -1.5  -4  -4   4 };
       * L_Minus_Minus is for the A1=-1, A2=-1 embedded adaptive intervention.;
	Eta_Plus_Plus = L_Plus_Plus*beta;
        * vector of fitted inverse-logit probabilities of Y=1 at each time point for (+1,+1);
    Eta_Plus_Minus = L_Plus_Minus*beta;
	    * vector of fitted inverse-logit probabilities of Y=1 at each time point for (+1,-1);
	Eta_Minus_Plus = L_Minus_Plus*beta;
        * vector of fitted inverse-logit probabilities of Y=1 at each time point for (-1,+1);
	Eta_Minus_Minus = L_Minus_Minus*beta;
	    * vector of fitted inverse-logit probabilities of Y=1 at each time point for (-1,-1);
	Mu_Plus_Plus = EXP(Eta_Plus_Plus)/(1+EXP(Eta_Plus_Plus));
        * vector of fitted probabilities of Y=1 at each time point for (+1,+1);
    Mu_Plus_Minus = EXP(Eta_Plus_Minus)/(1+EXP(Eta_Plus_Minus));
        * vector of fitted probabilities of Y=1 at each time point for (+1,-1);
	Mu_Minus_Plus = EXP(Eta_Minus_Plus)/(1+EXP(Eta_Minus_Plus));
        * vector of fitted probabilities of Y=1 at each time point for (-1,+1);
	Mu_Minus_Minus = EXP(Eta_Minus_Minus)/(1+EXP(Eta_Minus_Minus));
        * vector of fitted probabilities of Y=1 at each time point for (-1,-1);
    Estimated_Area_Plus_Plus = {0.5 1.0 1.0 1.0 1.0 0.5}*Mu_Plus_Plus/5;
        * Estimated area under the curve for the A1=+1, A2=+1 embedded adaptive intervention,
	      divided by the length of the time interval.;
	Estimated_Area_Plus_Minus = {0.5 1.0 1.0 1.0 1.0 0.5}*Mu_Plus_Minus/5;
        * Estimated area under the curve for the A1=+1, A2=-1 embedded adaptive intervention,
	      divided by the length of the time interval.;
	Estimated_Area_Minus_Plus = {0.5 1.0 1.0 1.0 1.0 0.5}*Mu_Minus_Plus/5;
        * Estimated area under the curve for the A1=-1, A2=+1 embedded adaptive intervention,
	      divided by the length of the time interval.;
	Estimated_Area_Minus_Minus = {0.5 1.0 1.0 1.0 1.0 0.5}*Mu_Minus_Minus/5;
        * Estimated area under the curve for the A1=-1, A2=-1 embedded adaptive intervention,
	      divided by the length of the time interval.;
	Estimate_PP_Versus_PM = Estimated_Area_Plus_Plus - Estimated_Area_Plus_Minus; 
        * The estimated difference equals the difference in the estimates;
    Estimate_PP_Versus_MP = Estimated_Area_Plus_Plus - Estimated_Area_Minus_Plus;
	Estimate_PP_Versus_MM = Estimated_Area_Plus_Plus - Estimated_Area_Minus_Minus;
	Estimate_PM_Versus_MP = Estimated_Area_Plus_Minus - Estimated_Area_Minus_Plus;
	Estimate_PM_Versus_MM = Estimated_Area_Plus_Minus - Estimated_Area_Minus_Minus;
	Estimate_MP_Versus_MM = Estimated_Area_Minus_Plus - Estimated_Area_Minus_Minus;
    Estimates =  Estimate_PP_Versus_PM  //
	             Estimate_PP_Versus_MP  //
	             Estimate_PP_Versus_MM  //
	             Estimate_PM_Versus_MP  // 
	             Estimate_PM_Versus_MM  //
	             Estimate_MP_Versus_MM ;
	* These are the choose(4,2) = 6 contrasts of interest among the AUC's for the
    * four embedded adaptive interventions.;
	Names = "AUC, ++ versus +-"  //
            "AUC, ++ versus -+"  //
            "AUC, ++ versus --"  //
            "AUC, +- versus -+"  //
            "AUC, +- versus --"  //
            "AUC, -+ versus --" ;
	Derivatives_PP_versus_PM = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Plus_Plus  // Mu_Plus_Minus)`#(1-(Mu_Plus_Plus  // Mu_Plus_Minus)`)/5;
	  * These are used in Cramer's delta method for converting the covariance in the 
	  * GEE coefficients to the variance of a contrast of interest.  # means 
	  * elementwise (one by one) multiplication.;
	Derivatives_PP_versus_MP = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Plus_Plus  // Mu_Minus_Plus)`#(1-(Mu_Plus_Plus  // Mu_Minus_Plus)`)/5;
	Derivatives_PP_versus_MM = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Plus_Plus  // Mu_Minus_Minus)`#(1-(Mu_Plus_Plus  // Mu_Minus_Minus)`)/5;
	Derivatives_PM_versus_MP = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Plus_Minus // Mu_Minus_Plus)`#(1-(Mu_Plus_Minus // Mu_Minus_Plus)`)/5;
	Derivatives_PM_versus_MM = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Plus_Minus // Mu_Minus_Minus)`#(1-(Mu_Plus_Minus // Mu_Minus_Minus)`)/5;
	Derivatives_MP_versus_MM = {.5 1 1 1 1 .5 -.5 -1 -1 -1 -1 -.5}#(Mu_Minus_Plus // Mu_Minus_Minus)`#(1-(Mu_Minus_Plus // Mu_Minus_Minus)`)/5; 
	Multipliers_PP_Versus_PM = L_Plus_Plus  // L_Plus_Minus;
	Multipliers_PP_Versus_MP = L_Plus_Plus  // L_Minus_Plus;
	Multipliers_PP_Versus_MM = L_Plus_Plus  // L_Minus_Minus;
	Multipliers_PM_Versus_MP = L_Plus_Minus // L_Minus_Plus;
	Multipliers_PM_Versus_MM = L_Plus_Minus // L_Minus_Minus;
	Multipliers_MP_Versus_MM = L_Minus_Plus // L_Minus_Minus; 
	Stderr_PP_Versus_PM = SQRT(  Derivatives_PP_Versus_PM  * 
	                              Multipliers_PP_Versus_PM * 
	                              CovBeta * 
	                              Multipliers_PP_Versus_PM` * 
	                              Derivatives_PP_Versus_PM`);
	   * This is an example of Cramer's delta method;
	Stderr_PP_Versus_MP = SQRT(  Derivatives_PP_Versus_MP * 
	                              Multipliers_PP_Versus_MP * 
	                              CovBeta * 
	                              Multipliers_PP_Versus_MP` * 
	                              Derivatives_PP_Versus_MP`);
	Stderr_PP_Versus_MM = SQRT(  Derivatives_PP_Versus_MM * 
	                              Multipliers_PP_Versus_MM * 
	                              CovBeta * 
	                              Multipliers_PP_Versus_MM` * 
	                              Derivatives_PP_Versus_MM`);
	Stderr_PM_Versus_MP = SQRT(  Derivatives_PM_Versus_MP * 
	                              Multipliers_PM_Versus_MP * 
	                              CovBeta * 
	                              Multipliers_PM_Versus_MP` * 
	                              Derivatives_PM_Versus_MP`);
	Stderr_PM_Versus_MM = SQRT(  Derivatives_PM_Versus_MM * 
	                              Multipliers_PM_Versus_MM * 
	                              CovBeta * 
	                              Multipliers_PM_Versus_MM` * 
	                              Derivatives_PM_Versus_MM`);
	Stderr_MP_Versus_MM = SQRT(  Derivatives_MP_Versus_MM * 
	                              Multipliers_MP_Versus_MM * 
	                              CovBeta * 
	                              Multipliers_MP_Versus_MM` * 
	                              Derivatives_MP_Versus_MM`);
 	StandardErrors = 	Stderr_PP_Versus_PM   //
				Stderr_PP_Versus_MP  // 
				Stderr_PP_Versus_MM  //
				Stderr_PM_Versus_MP  // 
				Stderr_PM_Versus_MM  // 
				Stderr_MP_Versus_MM  ;
	Z_Ratios = Estimates / StandardErrors;
    p_Values = 2*(1-PROBNORM(ABS(Z_Ratios)));
	CREATE Results VAR {Names Estimates StandardErrors Z_Ratios p_Values};
		APPEND ;
	CLOSE Results;
QUIT;

PROC PRINT DATA=Results; RUN;
