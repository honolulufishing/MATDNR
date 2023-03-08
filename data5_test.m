

function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt]=data5_test()
% Power Flow Data PJM 5-bus system Test System
%
%   Based on data from ...
%     F.Li and R.Bo, "Small Test Systems for Power System Economic Studies",
%     Proceedings of the 2010 IEEE Power & Energy Society General Meeting

%   Created by Rui Bo in 2006, modified in 2010, 2014.
%   Distributed with permission.

% Author(s):Chao Lei
% $Date:2023/02/02 10:35$

% IEEE 5-BUS TEST SYSTEM 
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang   ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg.  Min  Max
bus = [ 
    1	3	   0    0	    0	0		1	0	 0.90   1.1
	2	2	   0	    0	    0	0		1	0	 0.90   1.1
	3	1	   -500	    -98.61	    0	0		1	0	 0.90   1.1
	4	2	    0    0	    0	0		1	0	 0.90   1.1
	5	2	   0     0	    0	0		1	0	 0.90   1.1
];



    
% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol     gen     ---P---    gen_cost_b
%       NO.   MW   MVA  Max  Min  Mag.   status    max min     
gen = [ 
    1	0	    0	350	    -350	 1.0	    1  200 0  40
	2	100	    0	357.5	-357.5	 1.0	    1   200 0   20
	4	200	    0	390	-390	     1	        1   600 0  14
	5	100	    0	450	-450	     1.00	    1  600 0  20
      ];
% Note: sequence should be like ref+pv
% gen_cost_a = [2;2;2;2] ;   


% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R      X     B     switch   Cap.
%       NO.   from   to   p.u.   p.u.    p.u.     status   MW     
Lnbr_all = [
        1   1	5	0.00297	0.0297	0.00674  1  1000
        2   1	4	0.00297	0.0297	0.00674  1   1000
        3	2	3	0.00281	0.0281	0.00712  1  2000
        4   2	1	0.00304	0.0304	0.00658  1  5000
        5   5	2	0.00064	0.0064	0.03126 1  1000
        6   3	4	0.00108	0.0108	0.01852  1 100
];


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
         3 0    0.5   0      16 1 
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
bus (:,3:6) = bus (:,3:6)./sysdt(BASEMVA);
gen (:,2:5) = gen (:,2:5)./sysdt(BASEMVA);
Lnbr_all(:,8) = Lnbr_all(:,8)./sysdt(BASEMVA);
gen(:,8:9) = gen(:,8:9)./sysdt(BASEMVA);




