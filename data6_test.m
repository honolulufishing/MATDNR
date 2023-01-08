function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt]=data6_test()

% Power Flow Data IEEE 6-BUS Test System(Distributed Networks)
%
% Author(s):Chao Lei
% $Date:2013/03/01 10:35$


% IEEE 6-BUS TEST SYSTEM (American Electric Power)
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang   ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg.  Min  Max
bus = [   
       1      3    0        0    0   0    1    0     0.9   1.1
    2      1    -0.5       -0.05      0   0    1   0     0.9   1.1 
    3      1    -0.3    -0.03    0   0    1   0     0.9   1.1   
    4      1    -0.2    -0.02    0   0    1    0     0.9   1.1 
    5	   1	 -0.6	 -0.06	     0	 0	  1.150      0   0.90   1.1
     6	   1	 -0.1	  -0.01	     0	 0	  1.144     0.000     0.90   1.10
    ];

% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol
%       NO.   MW   MVA  Max  Min  Mag.
gen = [ 
	  	1	 0.5   0	1.5  -0.2	   1.1  1
      ];
% Note: sequence should be like ref+pv
  
  
% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R      X     1/2 B     switch
%       NO.   from   to   p.u.   p.u.    p.u.     status         
Lnbr_all = [
        1	   1	2	0.041	0.052	0   2
        2	   2	3	0.049	0.055	0   1
        3	   3	4	0.066	0.070	0   2
        4	   2	5	0.036	0.049	0   1
        5	   3	6	0.056	0.078	0   0   
        6      6    5   0.016	0.078	0   1
    ];
% 2 indicates this branch does not have a circuit breaker


% --------------------------- Transformer Data --------------------------- %
%      TrsBR.  Bus   Bus    R     X      ----- Tap -----     Tap      Tap
%       NO.   from   to    p.u.  p.u.    Ratio  Max  Min    Series   Status
trsfm = [       
 ]; 

% ---- Additonal Facility Data ---- %
% ----------- shunt capacity Data ------------ %
%      Bus  shtc     ---Q---  
%       NO. Mvar  Max  Min 
shtc = [   
         3 0.213    0.5   0      16 1 
]; 
% only reactive compensation exists in the PQ buses.

shtr= [];vctr=[];

% System Parameters
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
      ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

% Define Solver Parameters
sysdt(BASEMVA)   =   100;     % Base KV
sysdt(PFMAXIT)   =   20;     % Maximum Iteration
sysdt(ACCURACY)  =   1e-6;   % Tolerance Constant 
% sysdt(PFMETHOD)  =   2;    % Solver Method Index  
% sysdt(CPOPT)     =   1;    %  

% convert R and X in Line Branch with per unit system
%bus (:,3:6) = -bus (:,3:6);
%gen (:,2:3) = gen (:,2:3)./sysdt(BASEMVA);





