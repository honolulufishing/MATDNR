function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt]=data14()
% Power Flow Data IEEE 14-BUS Test System
%
% Author(s):Chao Lei
% $Date:2012/05/04 16:30$


% IEEE 14-BUS TEST SYSTEM (American Electric Power)
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang  ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg. Min  Max
bus = [ 1	   3	 0	     0	    0	0    1.06	0	0.9   1.1
	    2	   2	-21.7	-12.7	0	0	 1.045	0	0.9   1.1	
	    3	   2	-94.2	-19	    0	0	 1.01	0	0.9   1.1	
	    4	   1	-47.8	 3.9	0	0	 1.0 	0	0.9   1.1	
	    5	   1	-7.6	-1.6    0	0	 1.0 	0	0.9   1.1	
	    6	   2	-11.2	-7.5    0	0	 1.07	0	0.9   1.1	
	    7	   1	 0	     0	    0	0	 1.0 	0	0.9   1.1	
	    8	   2	 0	     0	    0	0	 1.09	0	0.9   1.1	
	    9	   1	-29.5	-16.6	0	0	 1.0	0	0.9   1.1	
	   10	   1	-9	    -5.8   	0	0	 1.0  	0	0.9   1.1	
	   11	   1	-3.5	-1.8   	0	0	 1.0 	0	0.9   1.1	
	   12	   1	-6.1    -1.6	0	0	 1.0 	0	0.9   1.1	
	   13	   1	0	0	0	0	 1.0	0		0.9   1.1
	   14	   1	-14.9	-5	    0	0	 1.0	0  0.9   1.1
       ];

% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol  gen     ---P---    gen_cost_b
%       NO.   MW   MVA  Max  Min  Mag. status    max min   
gen = [
       	1	 20	0	10	-20	   1.06	  1   575.88	0 10
	    2	40	0	50	-40	   1.045  1	100	0 17.5
	    3	 20	0	40	0	   1.01	1	140	0 18.5
	    6	 20	0	24	-6	   1.07	1	200	0 17
	    8	 20	0	24	-6	   1.09	1	350	0 11
      ];

% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R     X     1/2 B   switch   Cap.
%       NO.  from   to    p.u.  p.u.    p.u.      status   MW        
Lnbr_all = [
       1.   1	2	0.01938	0.05917	0.0528	1    1000
	   2.   1	5	0.05403	0.22304	0.0492	1    1000
	   3.   2	3	0.04699	0.19797	0.0438	1  1000
	   4.   2	4	0.05811	0.17632	0.034	1  1000
	   5.   2	5	0.05695	0.17388	0.0346	1  1000
	   6.   3	4	0.06701	0.17103	0.0128	1   1000
	   7.   4	5	0.01335	0.04211	0	1 1000
	   8.   6	11	0.09498	0.1989	0	1 1000
	   9.   6	12	0.12291	0.25581	0	1 1000
	  10.   6	13	0.06615	0.13027	0	1 1000
	  11.   9	10	0.03181	0.0845	0	1 1000
	  12.   9	14	0.12711	0.27038	0	1 1000
	  13.  10	11	0.08205	0.19207	0	1 1000
	  14.  12	13	0.22092	0.19988	0	1 1000
	  15.  13	14	0.17093	0.34802	0	1 1000
    ];

% --------------------------- Transformer Data --------------------------- %
%      TrsBR.  Bus   Bus    R     X      ----- Tap -----     Tap      Tap
%       NO.   from   to    p.u.  p.u.    Ratio  Max  Min    Series   Status
trsfm = [
         16.   4	7	0	0.20912 	 1	1.1  0.9      0         1
         17.   4	9	0	0.55618 	 1	1.1  0.9      0         1
         18.   5	6	0	0.25202 	 1	1.1  0.9      0         1
         19.   7	8	0	0.17615      1	    1.1  0.9      0         1
         20.   7	9	0	0.11001 	 1	    1.1  0.9      0         1
];

% ---- Additonal Facility Data ---- %
shtr= [];vctr=[];

% ----------- shunt capacity Data ------------ %
%      Bus  shtc     ---Q---  
%       NO. Mvar  Max  Min 
shtc = [   
         14  0.0   0.5   0      16 1 
]; 

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
