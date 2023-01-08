function [loss,gen,bus,Lnbr,trsfm,vctr,sysdt] = pfsoln(gen,bus,Lnbr,trsfm,vctr,sysdt,Ybus,Yf,Yt,V,ref,pv,pq,success,iterations)
% Newton-Raphson Power Flow with Standard Solution Structure  
%  [loss,gen,bus,Lnbr,trsfm,vctr,sysdt]=pfsoln(gen,bus,Lnbr,trsfm,vctr,sysdt,Ybus,Yf,Yt,V,ref,pv,pq,success,iterations) updates 
% bus, gen, branch data structures to match power flow soln.

% Author(s):Fang Liu, Chao Lei
% $Date:2012/05/04 22:32$

% constants
j = sqrt(-1);

% define named indices into bus, gen, branch matrices
[VCTR_BUS,VSET,WEIGHT,VCTR_VM,VCTR_OPTVM] = idx_VCTR;
[GEN_BUS,GEN_PG,GEN_QG,QGMAX,QGMIN,GEN_VG,GEN_OPTVG,GEN_OPTQG,GEN_OPTPG] = idx_Gen;
[BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus;
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
      ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;
[TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,TRS_PF,TRS_QF,TRS_PT,TRS_QT,...
      TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM;
[LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,LN_OPTPF,LN_OPTQF,...
      LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR;


% ----- update Qg for all gens and Pg for swing bus -----
% generator info
gbus = gen(:, GEN_BUS);					     % what buses are they at?
refgen = find(gen(:, GEN_BUS) == ref);		 % which is the reference gen?

% compute total injected bus powers
Sg = diag(V(gbus)) * conj(Ybus(gbus, :) * V);

% update Qg for all generators
gen(:, GEN_QG) = imag(Sg) - bus(gbus, QD);	 % inj Q + local Qd

% update Pg for swing bus
gen(refgen, GEN_PG) = real(Sg(refgen)) - bus(ref, PD);	 % inj P + local Pd

% the system total loss(active)
loss = sum(gen(:,GEN_PG)) + sum(bus(:,PD));                

% ----- update/compute branch power flows ----- %
if  ~isempty(trsfm)
    Sf = V([Lnbr(:, F_BUS);trsfm(:, F_BUS)]) .* conj(Yf * V);	% complex power at "from" bus
    St = V([Lnbr(:, T_BUS);trsfm(:, T_BUS)]) .* conj(Yt * V);	% complex power injected at "to" bus
else
    Sf = V([Lnbr(:, F_BUS)]) .* conj(Yf * V);	% complex power at "from" bus
    St = V([Lnbr(:, T_BUS)]) .* conj(Yt * V);	% complex power injected at "to" bus
end

% branch power flows 
nl = size(Lnbr, 1);		    % number of lines
nt = size(trsfm, 1);		% number of lines
j1 = nl + 1;              
j2 = nl + nt;
Lnbr(:, [LN_PF, LN_QF, LN_PT, LN_QT])      =   [ real(Sf(1:nl))  imag(Sf(1:nl))  real(St(1:nl))  imag(St(1:nl)) ];
trsfm(:,[TRS_PF, TRS_QF, TRS_PT, TRS_QT])  =   [ real(Sf(j1:j2)) imag(Sf(j1:j2)) real(St(j1:j2)) imag(St(j1:j2))];

Vm = abs(V);                              
Lnbr(:, [LN_FVM, LN_TVM])=[Vm(Lnbr(:,F_BUS)),Vm(Lnbr(:,T_BUS))];
trsfm(:,[TRS_FVM,TRS_TVM])=[Vm(trsfm(:,F_BUS)),Vm(trsfm(:,T_BUS))];
gen(:,GEN_VG) = Vm(gbus);

% ----- update bus voltages ----- %
bus(:, VM) = abs(V);
bus(:, VA) = angle(V);

if size(vctr,1) > 0  
    vctr(:,VCTR_VM) = bus(vctr(:,VCTR_BUS),VM);
end

sysdt(SUCCESS) = success;
sysdt(PFITER)  = iterations;


