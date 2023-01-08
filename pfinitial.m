function [Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr,trsfm,shtc,shtr,sysdt)
% Initialize Power Flow Analysis Package
% Note this program??package, we have considered all the shunts and capacity as constant loads.
% Author(s):Fang Liu, Chao Lei
% $Date:2012/05/04 17:30$

if nargin < 7
error('Please input 7 variables');
end
% -------------------- define named indices into bus, gen, branch matrices -------------------- %
% system indices
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;
% Generator indices
[GEN_BUS,GEN_PG,GEN_QG,QGMAX,QGMIN,GEN_VG,GEN_OPTVG,GEN_OPTQG,GEN_OPTPG] = idx_Gen;
% Bus indices
[BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus;
% Transformer indices
[TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,...
TRS_PF,TRS_QF,TRS_PT,TRS_QT,TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,...
TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM;
% Line Branch indices
[LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,...
LN_OPTPF,LN_OPTQF,LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR;
% Shunts indices
[SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC;
% Capacity indices
[SHTR_BUS,SHTR_QR,QRMAX,QRMIN,SHTR_SERIES,SHTR_STATE,OPTQR] = idx_SHTR;

% ------------- Pre-Calculation(including initial Voltage Mag, initial Sbus, Ybus) ------------- %
% System Parameters
pfmxt = sysdt(PFMAXIT); % Maximum iteration
tol = sysdt(ACCURACY); % Tolerance Accuracy
% Step 1: sort bus number according to bus types
nb = size(bus,1); % number of buses
Nbus = ones(nb,1);
% reference bus index
refgen = find(bus(:, BUS_TYPE) == 3);
ref = bus(refgen,BUS_I);
% PV bus indices
pvgen = find(bus(:, BUS_TYPE) == 2);
pv = bus(pvgen,BUS_I);
% PQ bus indices
Nbus([ref;pv]) = 0;
pq = find(Nbus);
% Step 2: Make Ybus
[Ybus,Yf,Yt] = makeYbus(bus,Lnbr,trsfm,shtc,shtr,pq,pv);
% Step 3: Calculate initial Sbus(all the negative initial load plus initial generation power)
gbus = gen(:,GEN_BUS);
Sbus = (bus(:,PD)+j*bus(:,QD)+bus(:,GS)+j*bus(:,BS));    % loads at all PQ buses
% loads plus generation at all PV and ref buses
Sbus(gbus) = Sbus(gbus)+gen(:,GEN_PG)+j*gen(:, GEN_QG);

if ~isempty(shtc)
shtcbus = shtc(:,SHTC_BUS);
Sbus(shtcbus) = Sbus(shtcbus) + j*shtc(:,SHTC_QC);
end

% Step 4: Initial bus magnitudes
V0 = ones(nb,1);
V0(gbus) = gen(:, GEN_VG);
% replace one with all the voltage magnitudes at PV and ref buses in the matrix V0


