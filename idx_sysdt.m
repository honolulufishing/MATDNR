function [DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
      ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt

DPRATE		= 1;    % depreciation rate per year 
INTRATE		= 2;	% interest rate per year 
OMRATE		= 3;	% operation maintenance cost rate
AUECOST		= 4;	% average unit electric cost 
AUCCOST		= 5;	% average unit capacitance cost (yuan/kvar)
AURCOST		= 6;	% average unit reactance cost (yuan/kvar)
TMAX	    = 7;	% max load consuming hour count (hours/year)
BASEMVA	    = 8;	% system base rating (MVA)
PFMETHOD	= 9;	% power flow method 
OPTMODEL	= 10;	% optimal model
OPTMETHOD	= 11;	% optimal method£¬=1£¬linear opt; =2, nonlinear opt; =3, genetic algorithm 
ACCURACY    = 12;	% termination tolerance
PFMAXIT		= 13;	% max iteration count
OPFMAXIT	= 14;	% max iteration count
POPNUM	    = 15;	% the number of population
CPOPT	    = 16;   % computer options  =1,power flow calculate; =2 var optimization computer
TARGET	    = 17;   % target function value
SUCCESS	    = 18;   % converged state
PFITER      = 19;   % the number of iteration in power flow
OPFITER     = 20;   % the number of iteration in reactive optimization 

return;
