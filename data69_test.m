function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data69_test()
% Power Flow Data IEEE 69-BUS Test System(Distributed Networks)
%
% Author(s):Chao Lei
% $Date:2020/03/14 22:35$


% IEEE 69-BUS TEST SYSTEM (American Electric Power)
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang   ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg.  Min  Max
bus = [   
1	3	0	0	0	0	1	0	0.9	1.1
2	1	30	3	0	0	1	0	0.9	1.1
3	1	30	3	0	0	1	0	0.9	1.1
4	1	30	3	0	0	1	0	0.9	1.1
5	1	13	0.3	0	0	1	0	0.9	1.1
6	1	26	2.2	0	0	1	0	0.9	1.1
7	1	40.4	30	0	0	1	0	0.9	1.1
8	1	75	54	0	0	1	0	0.9	1.1
9	1	30	22	0	0	1	0	0.9	1.1
10	1	28	19	0	0	1	0	0.9	1.1
11	1	145	104	0	0	1	0	0.9	1.1
12	1	145	104	0	0	1	0	0.9	1.1
13	1	8	5	0	0	1	0	0.9	1.1
14	1	8	5.5	0	0	1	0	0.9	1.1
15	1	30	3	0	0	1	0	0.9	1.1
16	1	45.5	30	0	0	1	0	0.9	1.1
17	1	60	35	0	0	1	0	0.9	1.1
18	1	60	35	0	0	1	0	0.9	1.1
19	1	30	3	0	0	1	0	0.9	1.1
20	1	10	0.6	0	0	1	0	0.9	1.1
21	1	114	81	0	0	1	0	0.9	1.1
22	1	50	3.5	0	0	1	0	0.9	1.1
23	1	30	0.3	0	0	1	0	0.9	1.1
24	1	28	20	0	0	1	0	0.9	1.1
25	1	30	0.3	0	0	1	0	0.9	1.1
26	1	14	10	0	0	1	0	0.9	1.1
27	1	14	10	0	0	1	0	0.9	1.1
28	1	26	18.6	0	0	1	0	0.9	1.1
29	1	26	18.6	0	0	1	0	0.9	1.1
30	1	30	3	0	0	1	0	0.9	1.1
31	1	30	3	0	0	1	0	0.9	1.1
32	1	30	3	0	0	1	0	0.9	1.1
33	1	14	10	0	0	1	0	0.9	1.1
34	1	9.5	14	0	0	1	0	0.9	1.1
35	1	60	4	0	0	1	0	0.9	1.1
36	1	26	18.55	0	0	1	0	0.9	1.1
37	1	26	18.55	0	0	1	0	0.9	1.1
38	1	30	0.3	0	0	1	0	0.9	1.1
39	1	24	17	0	0	1	0	0.9	1.1
40	1	24	17	0	0	1	0	0.9	1.1
41	1	12	1	0	0	1	0	0.9	1.1
42	1	30	3	0	0	1	0	0.9	1.1
43	1	60	4.3	0	0	1	0	0.9	1.1
44	1	30	3	0	0	1	0	0.9	1.1
45	1	39.22	26.3	0	0	1	0	0.9	1.1
46	1	39.22	26.3	0	0	1	0	0.9	1.1
47	1	30	0.3	0	0	1	0	0.9	1.1
48	1	79	56.4	0	0	1	0	0.9	1.1
49	1	384.7	274.5	0	0	1	0	0.9	1.1
50	1	384.7	274.5	0	0	1	0	0.9	1.1
51	1	40.5	28.3	0	0	1	0	0.9	1.1
52	1	36	2.7	0	0	1	0	0.9	1.1
53	1	4.35	3.5	0	0	1	0	0.9	1.1
54	1	26.4	19	0	0	1	0	0.9	1.1
55	1	24	17.2	0	0	1	0	0.9	1.1
56	1	30	3	0	0	1	0	0.9	1.1
57	1	30	3	0	0	1	0	0.9	1.1
58	1	30	3	0	0	1	0	0.9	1.1
59	1	100	72	0	0	1	0	0.9	1.1
60	1	30	3	0	0	1	0	0.9	1.1
61	1	1244	888	0	0	1	0	0.9	1.1
62	1	32	23	0	0	1	0	0.9	1.1
63	1	30	3	0	0	1	0	0.9	1.1
64	1	227	162	0	0	1	0	0.9	1.1
65	1	59	42	0	0	1	0	0.9	1.1
66	1	18	13	0	0	1	0	0.9	1.1
67	1	18	13	0	0	1	0	0.9	1.1
68	1	28	20	0	0	1	0	0.9	1.1
69	1	28	20	0	0	1	0	0.9	1.1

    ];

% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol
%       NO.   MW   MVA  Max  Min  Mag.
gen = [ 
	  	1	 0.5   0	10.5  -20	   1.1
      ];
% Note: sequence should be like ref+pv
  
  
% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R      X     1/2 B   switch
%       NO.   from   to   p.u.   p.u.    p.u.   status              
Lnbr_all = [
   1	1	2	0.0005	0.0012  0	2
   2	2	3	0.0005	0.0012	0   2
   3	3	4	0.0015	0.0036	0   1
   4	4	5	0.0251	0.0294  0	1
   5	5	6	0.3660	0.1864	0   1
   6	6	7	0.3811	0.1941	0   1
   7	7	8	0.0922	0.0470	0   1
   8	8	9	0.0493	0.0251	0   1
   9	9	10	0.8190	0.2707	0   1
   10	10	11	0.1872	0.0619	0   1
   11	11	12	0.7114	0.2351	0    1
12	12	13	1.0300	0.3400	0  0
13	13	14	1.0440	0.3450 0	1
14	14	15	1.0580	0.3496	0 1
15	15	16	0.1966	0.0650	0 1
16	16	17	0.3744	0.1238	0 1
17	17	18	0.0047	0.0016	0 1
18	18	19	0.3276	0.1083	0 1
19	19	20	0.2106	0.0690	0 1
20	20	21	0.3416	0.1129	0 1
21	21	22	0.0140	0.0046	0 1
22	22	23	0.1591	0.0526	0 1
23	23	24	0.3463	0.1145	0 1
24	24	25	0.7488	0.2475	0 1
25	25	26	0.3089	0.1021	0 1
26	26	27	0.1732	0.0572	0 1

27	3	28	0.0044	0.0108	0 2
28	28	29	0.0640	0.1565	0 2
29	29	30	0.3978	0.1315	0 2
30	30	31	0.0702	0.0232	0 2
31	31	32	0.3510	0.1160	0 2
32	32	33	0.8390	0.2816	0 2
33	33	34	1.7080	0.5646	0 2
34	34	35	1.4740	0.4873	0 2

35	3	36	0.0044	0.0108	0 1
36	36	37	0.0640	0.1565	0 1
37	37	38	0.1053	0.1230	0 1
38	38	39	0.0304	0.0355	0 1
39	39	40	0.0018	0.0021	0 1
40	40	41	0.7283	0.8509	0 1
41	41	42	0.3100	0.3623	0 1
42	42	43	0.0410	0.0478	0 1
43	43	44	0.0092	0.0116	0 1
44	44	45	0.1089	0.1373	0 1
45	45	46	0.0009	0.0012	0 1

46	4	47	0.0034	0.0084	0 1
47	47	48	0.0851	0.2083	0 1
48	48	49	0.2898	0.7091	0 1
49	49	50	0.0822	0.2011	0 1

50	8	51	0.0928	0.0473	0 2
51	51	52	0.3319	0.1114	0 2

52	9	53	0.1740	0.0886	0 1
53	53	54	0.2030	0.1034	0 1
54	54	55	0.2842	0.1447	0 1
55	55	56	0.2813	0.1433	0 1
56	56	57	1.5900	0.5337	0 1
57	57	58	0.7837	0.2630	0 0
58	58	59	0.3042	0.1006	0 1
59	59	60	0.3861	0.1172	0 1
60	60	61	0.5075	0.2585	0 1
61	61	62	0.0974	0.0496	0 1
62	62	63	0.1450	0.0738	0 0
63	63	64	0.7105	0.3619	0 1
64	64	65	1.0410	0.5302	0 1

65	11	66	0.2012	0.0611	0 2
66	66	67	0.0047	0.0014	0 2

67	12	68	0.7394	0.2444	0 2
68	68	69	0.0047	0.0016	0 2

69   27	65	0.0047	0.0016	0 1
70   50	59	0.0047	0.0016	0 1
71   15	46	0.0047	0.0016	0 1
72   11	43	0.0047	0.0016	0 0
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
    10 0 0 0
%     10 0.0 0.1 0
%     27 0.0 0.1 0
%     28 0.0 0.1 0
%     35 0.0 0.1 0
%     40 0.0 0.1 0
%     46 0.0 0.1 0
%     50 0.0 0.1 0
%     52 0.0 0.1 0
%     65 0.0 0.1 0
%     67 0.0 0.1 0
];
% only reactive compensation exists in the PQ buses.

shtr= [];vctr=[];

% System Parameters
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
      ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

% Define Solver Parameters
sysdt(BASEMVA)   =   1000;     % Base MVA, V_Base(kV) 12.66 Z_base(ohms) 1.60275
sysdt(PFMAXIT)   =   20;     % Maximum Iteration
sysdt(ACCURACY)  =   1e-6;   % Tolerance Constant 
% sysdt(PFMETHOD)  =   2;    % Solver Method Index  
% sysdt(CPOPT)     =   1;    %  

% convert R and X in Line Branch with per unit system
bus (:,3:6) = -bus (:,3:6)./sysdt(BASEMVA);
Lnbr_all (:,4:5) = Lnbr_all(:,4:5)./12.66^2;




