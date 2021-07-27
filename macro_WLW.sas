*************************************************************************************************************
*************************************************************************************************************
*************************************************************************************************************

WLW
WLW with optimal weights or combined Z-scores for Multiple Types of Events

Inputs:

Data		Name of the dataset.

ID			Name of the variable for identifying the subject so that all observations of the same subject 
			have the same ID value.

Enum        Name of the variable to index the multiple events. For example, Enum=1 for the first event, 
			Enum=2 for the second event, and so on. Indices can be either numbers or characters.
			
Time		Name of the variable to represent the observed time from some time origin for the event.

Status		Name of the variable to indicate whether the Time value is a censored or uncensored time, 
			with the value that indicates censoring enclosed in parenthesis. Same as the MODEL statement 
			in SAS PHREG, Status=1 is the default for censoring and no parenthesis is needed. Otherwise, 
			for example, if Status=1 indicates an uncensored time and Status=0 indicates a censored time, 
			then Status(0) should be used.   
			
Treatment	Name of the binary treatment indicator. Treatment variable should be coded in such a way that 
			a positive value of the coefficient means that the treatment is beneficial. 

Covariates	Names of independent prognostic factors. Must be space delimited if there are more than one 
			factors (e.g., Covariates=Var1 Var2 Var3). Prognostic factors must be numeric.

Test		Alternative hypothesis. Must be "one.sided" (default) or "two.sided". If Test=one.sided, a 
			one-sided test will be performed. If Test=two.sided, a two-sided test will be performed.	

*************************************************************************************************************
*************************************************************************************************************
************************************************************************************************************;

%macro WLW(       Data = , 
				    ID = , 
			      Enum = , 
			      Time = , 
			    Status = , 
		     Treatment = , 
		    Covariates = %str(NONE), 
		          Test = %str(one.sided)
		  );
	
	%global numtype numcov;

****************;
****************;
****************;

	proc sort data=%str(&Data);
		by &ID &Enum;
	run;

	proc freq data=%str(&Data) nlevels;
		ods exclude all;
		ods output NLevels=_work_numtype;
		tables &Enum;
	run;

	ods select all;

	data _null_;
		set _work_numtype;
		call symput("numtype", left(nlevels));
		call symput("numcov", left(countw(compbl("&Covariates"), " ")));
	run;

	proc sql noprint;
		select distinct(&Enum) into :types separated by ' ' from %str(&Data);
		%let type = &types;
	quit;

	data _null_;
		%let newtrt = ;
		%do j=1 %to &numtype;
			%let newtrt = &newtrt &Treatment.%scan(&type, &j);
		%end;
		
		%let newcovars = ;
		%do i=1 %to &numcov;
			%let covar = %scan(&Covariates, &i);
			%do j=1 %to &numtype;
				%let newcovars = &newcovars &covar.%scan(&type, &j);
			%end;
		%end;
	run;
	
****************;
****************;
****************;

	proc phreg data=%str(&Data) covs(aggregate);
		model &Time*&Status=&newtrt &newcovars/covb;
		strata &Enum;
		id &ID;
		ods output CovB = _work_cov 
		           ParameterEstimates=_work_beta 
			       TestAverage=_work_average
			       TestCoeff = _work_testcoef
			       ParameterEstimates = _work_mle;

		%do j=1 %to &numtype;
			%let temptrt = %scan(&newtrt, &j);
			%let temptype = %scan(&type, &j);
			&temptrt=&Treatment * (&Enum="&temptype");
		%end;

		%do i=1 %to &numcov;
			%let covar = %scan(&Covariates, &i);

			%do j=1 %to &numtype;
				%let tempcovar = &covar.%scan(&type, &j);
				%let temptype = %scan(&type, &j);
				&tempcovar=&covar * (&Enum="&temptype");
			%end;
		%end;
		%let testtrt = %sysfunc(tranwrd(%quote(&newtrt), %str( ), %str(,)));
        Treatment: test &testtrt / e average;
	run;

****************;
****************;
****************;

*** Summary Outputs ***;

	data _work_mle;
		set _work_mle;
		drop StdErrRatio HazardRatio;
	run;
	
	data _work_testcoef;
		set _work_testcoef;
		if _n_ <=2;
		keep Parameter Average;
		label Average="Optimal weight";
	run;

	proc iml;
		mle = TableCreateFromDataSet("_work_mle");
		call TablePrint(mle) label="Analysis of Maximum Likelihood Estimates" 
		          			 COLHEADER="Names"
			                 ID="Parameter";
		
		coef = TableCreateFromDataSet("_work_testcoef");
		call TablePrint(coef) label="Optimal Weights for Test Treatment"
		                      COLHEADER="Labels"
		                      ID="Parameter";
		
		use _work_beta;
		read all var{Estimate} into beta;
		close;
		use _work_cov;
		read all into cov;
		close;
		beta = beta[1:2];
		cov = cov[1:2,1:2];
		se=sqrt(vecdiag(cov));
		w=1/se;
		T_nu=t(w)*beta;
		T_de=sqrt(t(w)*cov*w);
		T=T_nu/T_de;

		%if %lowcase(%qsysfunc(dequote(&Test)))=one.sided %then %do;
				use _work_average;
				read all var{zScore} into T_ave;
				close;
				pval=1-cdf('normal', T);
				pval_ave=1-cdf('normal', T_ave);
				method={"Optimal weights", "Combined Z-scores"};
				T_all={0, 0};
				T_all[1]=T_ave;
				T_all[2]=T;
				P_all={0, 0};
				P_all[1]=pval_ave;
				P_all[2]=pval;
				tbl=TableCreate("Method", method);
				call TableAddVar(tbl, {"Test statistic" "P-value"}, T_all||P_all);
				call TablePrint(tbl) 
					label="One-Sided Test Results of the WLW Method" var={"Method" 
					"Test statistic" "P-value"} justify={'L' 'C' 'C'} ID="Method";
		%end;

		%if %lowcase(%qsysfunc(dequote(&Test)))=two.sided %then %do;
				use _work_average;
				read all var{zScore} into T_ave;
				read all var{Probz} into pval_ave;
				close;
				pval=2*(1-cdf('normal', abs(T)));
				method={"Optimal weights", "Combined Z-scores"};
				T_all={0, 0};
				T_all[1]=T_ave;
				T_all[2]=T;
				P_all={0, 0};
				P_all[1]=pval_ave;
				P_all[2]=pval;
				tbl=TableCreate("Method", method);
				call TableAddVar(tbl, {"Test statistic" "P-value"}, T_all||P_all);
				call TablePrint(tbl) 
					label="Two-Sided Test Results of the WLW Method" var={"Method" 
					"Test statistic" "P-value"} justify={'L' 'C' 'C'} ID="Method";
		%end;
	quit;
	
%mend WLW;

************************************************************************************************************;
************************************************************************************************************;
************************************************************************************************************;

data example1;
	infile 'example1.dat';
	input id treatment basecat event time status;
	basecat4 = (basecat=4);
	basecat5 = (basecat=5);
run;

data example2;
	infile 'example2.dat';
	input id treatment basecat event time status;
	basecat4 = (basecat=4);
	basecat5 = (basecat=5);
run;

%WLW(Data=example1, ID=id, Enum=event, Time=time, Status=status(0), Treatment=treatment, 
	Covariates=basecat4 basecat5, Test=one.sided)
%WLW(Data=example2, ID=id, Enum=event, Time=time, Status=status(0), Treatment=treatment, 
	Covariates=basecat4 basecat5, Test=one.sided)
