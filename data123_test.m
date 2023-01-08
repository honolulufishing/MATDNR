function [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data123_test()
% Power Flow Data IEEE 123-BUS Test System (Distributed Networks)
%
% Author(s):Chao Lei
% $Date:2022/08/05 22:35$


% IEEE 123-BUS TEST SYSTEM
%
% -------------------- Bus Data ----------------- %
%        Bus  Bus   ---Load---    -Injected-  Vol   Ang   ---Vol---
%        No   Type  MW    Mvar     GS   Bs    Mag.  Deg.  Min  Max
bus = [   
1	3	0	0	0	0	1	0	0.9	1.1
2	1	20	10	0	0	1	0	0.9	1.1
3	1	0	0	0	0	1	0	0.9	1.1
4	1	40	20	0	0	1	0	0.9	1.1
5	1	20	10	0	0	1	0	0.9	1.1
6	1	40	20	0	0	1	0	0.9	1.1
7	1	20	10	0	0	1	0	0.9	1.1
8	1	0	0	0	0	1	0	0.9	1.1
9	1	40	20	0	0	1	0	0.9	1.1
10	1	20	10	0	0	1	0	0.9	1.1
11	1	40	20	0	0	1	0	0.9	1.1
12	1	20	10	0	0	1	0	0.9	1.1
13	1	0	0	0	0	1	0	0.9	1.1
14	1	0	0	0	0	1	0	0.9	1.1
15	1	0	0	0	0	1	0	0.9	1.1
16	1	40	20	0	0	1	0	0.9	1.1
17	1	20	10	0	0	1	0	0.9	1.1
18	1	0	0	0	0	1	0	0.9	1.1
19	1	40	20	0	0	1	0	0.9	1.1
20	1	40	20	0	0	1	0	0.9	1.1
21	1	0	0	0	0	1	0	0.9	1.1
22	1	40	20	0	0	1	0	0.9	1.1
23	1	0	0	0	0	1	0	0.9	1.1
24	1	40	20	0	0	1	0	0.9	1.1
25	1	0	0	0	0	1	0	0.9	1.1
26	1	0	0	0	0	1	0	0.9	1.1
27	1	0	0	0	0	1	0	0.9	1.1
28	1	40	20	0	0	1	0	0.9	1.1
29	1	40	20	0	0	1	0	0.9	1.1
30	1	40	20	0	0	1	0	0.9	1.1
31	1	20	10	0	0	1	0	0.9	1.1
32	1	20	10	0	0	1	0	0.9	1.1
33	1	40	20	0	0	1	0	0.9	1.1
34	1	40	20	0	0	1	0	0.9	1.1
35	1	40	20	0	0	1	0	0.9	1.1
36	1	0	0	0	0	1	0	0.9	1.1
37	1	40	20	0	0	1	0	0.9	1.1
38	1	20	10	0	0	1	0	0.9	1.1
39	1	20	10	0	0	1	0	0.9	1.1
40	1	0	0	0	0	1	0	0.9	1.1
41	1	20	10	0	0	1	0	0.9	1.1
42	1	20	10	0	0	1	0	0.9	1.1
43	1	40	20	0	0	1	0	0.9	1.1
44	1	0	0	0	0	1	0	0.9	1.1
45	1	20	10	0	0	1	0	0.9	1.1
46	1	20	10	0	0	1	0	0.9	1.1
47	1	105	75	0	0	1	0	0.9	1.1
48	1	210	150	0	0	1	0	0.9	1.1
49	1	140	95	0	0	1	0	0.9	1.1
50	1	40	20	0	0	1	0	0.9	1.1
51	1	20	10	0	0	1	0	0.9	1.1
52	1	40	20	0	0	1	0	0.9	1.1
53	1	40	20	0	0	1	0	0.9	1.1
54	1	0	0	0	0	1	0	0.9	1.1
55	1	20	10	0	0	1	0	0.9	1.1
56	1	20	10	0	0	1	0	0.9	1.1
57	1	0	0	0	0	1	0	0.9	1.1
58	1	20	10	0	0	1	0	0.9	1.1
59	1	20	10	0	0	1	0	0.9	1.1
60	1	20	10	0	0	1	0	0.9	1.1
61	1	0	0	0	0	1	0	0.9	1.1
62	1	40	20	0	0	1	0	0.9	1.1
63	1	40	20	0	0	1	0	0.9	1.1
64	1	75	35	0	0	1	0	0.9	1.1
65	1	140	100	0	0	1	0	0.9	1.1
66	1	75	35	0	0	1	0	0.9	1.1
67	1	0	0	0	0	1	0	0.9	1.1
68	1	20	10	0	0	1	0	0.9	1.1
69	1	40	20	0	0	1	0	0.9	1.1
70	1	20	10	0	0	1	0	0.9	1.1
71	1	40	20	0	0	1	0	0.9	1.1
72	1	0	0	0	0	1	0	0.9	1.1
73	1	40	20	0	0	1	0	0.9	1.1
74	1	40	20	0	0	1	0	0.9	1.1
75	1	40	20	0	0	1	0	0.9	1.1
76	1	245	180	0	0	1	0	0.9	1.1
77	1	40	20	0	0	1	0	0.9	1.1
78	1	0	0	0	0	1	0	0.9	1.1
79	1	40	20	0	0	1	0	0.9	1.1
80	1	40	20	0	0	1	0	0.9	1.1
81	1	0	0	0	0	1	0	0.9	1.1
82	1	40	20	0	0	1	0	0.9	1.1
83	1	20	10	0	0	1	0	0.9	1.1
84	1	20	10	0	0	1	0	0.9	1.1
85	1	40	20	0	0	1	0	0.9	1.1
86	1	20	10	0	0	1	0	0.9	1.1
87	1	40	20	0	0	1	0	0.9	1.1
88	1	40	20	0	0	1	0	0.9	1.1
89	1	0	0	0	0	1	0	0.9	1.1
90	1	40	20	0	0	1	0	0.9	1.1
91	1	0	0	0	0	1	0	0.9	1.1
92	1	40	20	0	0	1	0	0.9	1.1
93	1	0	0	0	0	1	0	0.9	1.1
94	1	40	20	0	0	1	0	0.9	1.1
95	1	20	10	0	0	1	0	0.9	1.1
96	1	20	10	0	0	1	0	0.9	1.1
97	1	0	0	0	0	1	0	0.9	1.1
98	1	40	20	0	0	1	0	0.9	1.1
99	1	40	20	0	0	1	0	0.9	1.1
100	1	40	20	0	0	1	0	0.9	1.1
101	1	0	0	0	0	1	0	0.9	1.1
102	1	20	10	0	0	1	0	0.9	1.1
103	1	40	20	0	0	1	0	0.9	1.1
104	1	40	20	0	0	1	0	0.9	1.1
105	1	0	0	0	0	1	0	0.9	1.1
106	1	40	20	0	0	1	0	0.9	1.1
107	1	40	20	0	0	1	0	0.9	1.1
108	1	0	0	0	0	1	0	0.9	1.1
109	1	40	20	0	0	1	0	0.9	1.1
110	1	0	0	0	0	1	0	0.9	1.1
111	1	20	10	0	0	1	0	0.9	1.1
112	1	20	10	0	0	1	0	0.9	1.1
113	1	40	20	0	0	1	0	0.9	1.1
114	1	20	10	0	0	1	0	0.9	1.1
115	1	0	0	0	0	1	0	0.9	1.1
116	1	0	0	0	0	1	0	0.9	1.1
117	1	0	0	0	0	1	0	0.9	1.1
118	1	0	0	0	0	1	0	0.9	1.1
119	1	0	0	0	0	1	0	0.9	1.1
120	1	0	0	0	0	1	0	0.9	1.1
121	1	0	0	0	0	1	0	0.9	1.1
122	1	0	0	0	0	1	0	0.9	1.1
123	1	40	20	0	0	1	0	0.9	1.1
    ];

% -------------------- Gen Data ----------------- %
%      Bus    --Gen--   ---Q---   Vol
%       NO.   MW   MVA  Max  Min  Mag.
gen = [ 
	  	1	 0   0	10.5  -20	   1.1  1
      ];
% Note: sequence should be like ref+pv
  
  
% -------------------- Line Data ----------------- %
%      LnBR.  Bus   Bus    R      X     1/2 B   switch
%       NO.   from   to   p.u.   p.u.    p.u.   status              
Lnbr_all = [
1	1	123	0.002	0.0047	0	2
2	123	2	0.0025	0.0026	0	2
3	123	3	0.0036	0.0037	0	2
4	123	7	0.0015	0.0035	0	2
5	3	4	0.0029	0.0029	0	2
6	3	5	0.0047	0.0048	0	2
7	5	6	0.0036	0.0037	0	2
8	7	8	0.001	0.0024	0	2
9	8	12	0.0033	0.0033	0	2
10	8	9	0.0033	0.0033	0	2
11	8	13	0.0015	0.0035	0	2
12	9	14	0.0062	0.0063	0	2
13	13	34	0.0022	0.0022	0	2
14	14	11	0.0036	0.0037	0	2
15	14	10	0.0036	0.0037	0	2
16	15	16	0.0055	0.0055	0	2
17	15	17	0.0051	0.0052	0	2
18	18	19	0.0036	0.0037	0	2
19	18	21	0.0015	0.0034	0	2
20	19	20	0.0047	0.0048	0	2
21	21	22	0.0076	0.0077	0	2
22	21	23	0.0013	0.0029	0	2
23	23	24	0.008	0.0081	0	2
24	23	25	0.0014	0.0032	0	2
25	25	26	0.0018	0.0041	0	2
26	25	28	0.001	0.0023	0	2
27	26	27	0.0014	0.0032	0	2
28	26	31	0.0033	0.0033	0	2
29	27	33	0.0073	0.0074	0	2
30	28	29	0.0015	0.0034	0	2
31	29	30	0.0018	0.004	0	2
32	30	122	0.001	0.0023	0	2
33	31	32	0.0044	0.0044	0	2
34	34	15	0.0015	0.0015	0	2
35	35	36	0.0033	0.0077	0	2
36	36	37	0.0044	0.0044	0	2
37	36	38	0.0036	0.0037	0	2
38	38	39	0.0047	0.0048	0	2
39	40	41	0.0047	0.0048	0	2
40	42	43	0.0073	0.0074	0	2
41	44	45	0.0029	0.0029	0	2
42	45	46	0.0044	0.0044	0	2
43	47	48	0.0008	0.0017	0	2
44	54	55	0.0014	0.0032	0	2
45	55	56	0.0014	0.0032	0	2
46	57	58	0.0036	0.0037	0	2
47	58	59	0.0036	0.0037	0	2
48	60	61	0.0028	0.0063	0	2
49	60	62	0.0042	0.0021	0	2
50	62	63	0.0029	0.0014	0	2
51	63	64	0.0058	0.0029	0	2
52	64	65	0.0071	0.0035	0	2
53	65	66	0.0054	0.0027	0	2
54	67	68	0.0029	0.0029	0	2
55	68	69	0.004	0.0041	0	2
56	69	70	0.0047	0.0048	0	2
57	70	71	0.004	0.0041	0	2
58	72	73	0.004	0.0041	0	2
59	73	74	0.0051	0.0052	0	2
60	74	75	0.0058	0.0059	0	2
61	76	77	0.002	0.0047	0	2
62	77	78	0.0005	0.0012	0	2
63	78	79	0.0011	0.0027	0	2
64	78	80	0.0024	0.0056	0	2
65	80	81	0.0024	0.0056	0	2
66	81	82	0.0013	0.0029	0	2
67	81	84	0.0098	0.01	0	2
68	82	83	0.0013	0.0029	0	2
69	84	85	0.0069	0.007	0	2
70	87	88	0.0025	0.0026	0	2
71	89	90	0.0033	0.0033	0	2
72	91	92	0.0044	0.0044	0	2
73	93	95	0.0015	0.0035	0	2
74	95	96	0.0029	0.0029	0	2
75	97	98	0.0014	0.0032	0	2
76	98	99	0.0028	0.0064	0	2
77	99	100	0.0015	0.0035	0	2
78	100	118	0.004	0.0093	0	2
79	101	102	0.0033	0.0033	0	2
80	102	103	0.0047	0.0048	0	2
81	103	104	0.0102	0.0103	0	2
82	105	106	0.0033	0.0033	0	2
83	106	107	0.0084	0.0085	0	2
84	108	109	0.0065	0.0066	0	2
85	109	110	0.0044	0.0044	0	2
86	110	111	0.0084	0.0085	0	2
87	110	112	0.0018	0.0018	0	2
88	112	113	0.0076	0.0077	0	2
89	113	114	0.0047	0.0048	0	2
90	13	18	0.0042	0.0095	0	1
91	18	121	0.0013	0.0029	0	1
92	121	35	0.0019	0.0044	0	1
93	35	40	0.0013	0.0029	0	1
94	40	42	0.0013	0.0029	0	1
95	42	44	0.001	0.0024	0	1
96	44	47	0.0013	0.0029	0	1
97	47	49	0.0013	0.0029	0	1
98	49	50	0.0013	0.0029	0	1
99	50	51	0.0013	0.0029	0	1
100	51	116	0.0025	0.0058	0	1
101	116	115	0.0013	0.0029	0	1
102	115	108	0.0051	0.0117	0	1
103	108	105	0.0016	0.0038	0	1
104	105	101	0.0014	0.0032	0	0
105	101	117	0.0013	0.0029	0	1
106	117	97	0.0013	0.0029	0	1
107	97	67	0.0013	0.0029	0	1
108	13	120	0.0013	0.0029	0	1
109	120	52	0.002	0.0047	0	1
110	52	53	0.001	0.0024	0	1
111	53	54	0.0006	0.0015	0	1
112	54	57	0.0018	0.0041	0	1
113	57	60	0.0038	0.0087	0	1
114	60	119	0.0013	0.0029	0	1
115	119	67	0.0018	0.0041	0	1
116	54	94	0.0013	0.0029	0	1
117	94	93	0.004	0.0041	0	1
118	93	91	0.0011	0.0027	0	1
119	91	89	0.0011	0.0027	0	1
120	89	87	0.0014	0.0032	0	1
121	87	86	0.0023	0.0053	0	1
122	86	76	0.0035	0.0082	0	1
123	76	72	0.001	0.0023	0	0
124	72	67	0.0014	0.0032	0	1
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
12	 0.1667 0.9	-0.3 16 1 
35	 0.5633 0.9	-0.3 16 1 
54	 0.7808  0.9	-0.3 16 1 
76	0.3683 0.9	-0.3 16 1 
108 0.2324	0.9	-0.3 16 1 
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




