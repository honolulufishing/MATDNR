function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt]=data16_test()
% Power Flow Data IEEE 16-BUS Test System(Distributed Networks)
%
% Author(s):Chao Lei
% $Date:2022/08/14 10:35$

%   Data from ...
%       Civanlar S, Grainger JJ, Yin H, Lee SSH (1988) Distribution Feeder
%       Reconfiguration for Loss Reduction. IEEE Trans Power Deliv
%       3:1217-1223. doi: 10.1109/61.193906
%       URL: https://doi.org/10.1109/61.193906
%   and
%       Zhu JZ (2002) Optimal reconfiguration of electrical distribution
%       network using the refined genetic algorithm, Electr Power Syst Res
%       62:37-42. doi: 10.1016/S0378-7796(02)00041-X
%       URL: https://doi.org/10.1016/S0378-7796(02)00041-X


% IEEE 16-BUS TEST SYSTEM (American Electric Power)
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang   ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg.  Min  Max
bus = [   
    1	3	0	     0	        0	0	1	0	0.9	1.1;
	2	1	200	     100	        0	0	1	0	0.9	1.1;
    3	1	200	     100	        0	0	1	0	0.9	1.1;
	4	1	200	     100	        0	0	1	0	0.9	1.1;
	5	1	2000	1600		0	0	1	0	0.9	1.1;
	6	1	3000	400	        0	0	1	0	0.9	1.1;
	7	1	2000	400		0	0	1	0	0.9	1.1;
	8	1	1500	1200		0	0	1	0	0.9	1.1;
	9	1	4000	2700		0	0	1	0	0.9	1.1;
	10	1	5000	1800		0	0	1	0	0.9	1.1;
	11	1	1000	900	        0	0	1	0	0.9	1.1;
	12	1	600  	500	    0	0	1	0	0.9	1.1;
	13	1	4500	1700		0	0	1	0	0.9	1.1;
	14	1	1000	900	        0	0	1	0	0.9	1.1;
	15	1	1000	1100		0	0	1	0	0.9	1.1;
	16	1	1000	900	        0	0	1	0	0.9	1.1;
	17	1	2100	800		0	0	1	0	0.9	1.1;   
    ];

% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol
%       NO.   MW   MVA  Max  Min  Mag.
gen = [ 
	  	1	 0.5   0	2.5  -0.2	   1.1 1
      ];
% Note: sequence should be like ref+pv
  
  
% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R      X     1/2 B     switch
%       NO.   from   to   p.u.   p.u.    p.u.     status         
Lnbr_all = [
    1 1 2   0.075	0.1	        0	2 
    2 2	5	0.075	0.1	        0	1
	3 5	6	0.08	0.11	    0	1
	4 5	7	0.09	0.18	    0	1
	5 7	8	0.04	0.04	    0	1
	6 1 3  0.075	0.1	        0	1 
    7 3	9	0.11	0.11	    0	1
	8 9	10	0.08	0.11	    0	1
	9 9	11	0.11	0.11	    0	0
	10 10	12	0.11	0.11	0	0
	11 10	13	0.08	0.11	0	2
	12 1 4  0.075	0.1	        0	1 
    13 4	14	0.11	0.11	0	1
	14 14	15	0.09	0.12	0	1
	15 14	16	0.08	0.11	0	1
	16 16	17	0.04	0.04	0	1
	17 6	12	0.04	0.04	0	1
	18 11	15	0.04	0.04	0	1
	19 8	17	0.09	0.12	0	0
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
shtc = [ 13   0.473 0.6 0 16 1 
]; 
% only reactive compensation exists in the PQ buses.

shtr= [];vctr=[];

% System Parameters
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
      ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

% Define Solver Parameters
sysdt(BASEMVA)   =   10000;     % Base KV
sysdt(PFMAXIT)   =   20;     % Maximum Iteration
sysdt(ACCURACY)  =   1e-6;   % Tolerance Constant 
% sysdt(PFMETHOD)  =   2;    % Solver Method Index  
% sysdt(CPOPT)     =   1;    %  

% convert R and X in Line Branch with per unit system
bus (:,3:6) = -bus (:,3:6)./sysdt(BASEMVA);

Lnbr_all(:,4:5) = Lnbr_all(:,4:5)./(12.66^2/10);




