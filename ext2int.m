function [i2e,gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr)
% EXT2INT   Converts external to internal bus numbering.
%   converts external bus numbers (possibly non-consecutive) to consecutive internal
%   bus numbers which start at 1.

% --------------- define named indices into bus, gen, branch matrices --------------- %
[GEN_BUS,GEN_PG,GEN_QG,QGMAX,QGMIN,GEN_VG,GEN_OPTVG,GEN_OPTQG,GEN_OPTPG] = idx_Gen;
[BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus;
[TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,TRS_PF,TRS_QF,TRS_PT,TRS_QT,...
      TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM;
[LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,LN_OPTPF,LN_OPTQF,...
      LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR;
[SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC;
[SHTR_BUS,SHTR_QR,QRMAX,QRMIN,SHTR_SERIES,SHTR_STATE,OPTQR] = idx_SHTR;
[VCTR_BUS,VSET,WEIGHT,VCTR_VM,VCTR_OPTVM] = idx_VCTR;


i2e	= bus(:, BUS_I);
e2i = zeros(max(i2e), 1);
e2i(i2e) = [1:size(bus, 1)]'; 

bus(:, BUS_I)	= e2i( bus(:, BUS_I)   );
gen(:, GEN_BUS)	= e2i( gen(:, GEN_BUS) );

if size(shtc,1)>0
    shtc(:,SHTC_BUS) = e2i( shtc(:, SHTC_BUS) );
end
if size(shtr,1)>0
    shtr(:,SHTR_BUS) = e2i( shtr(:, SHTR_BUS) );
end
if size(vctr,1)>0
    vctr(:,VCTR_BUS)= e2i( vctr(:, VCTR_BUS) );
end

Lnbr(:, F_BUS)	= e2i( Lnbr(:, F_BUS) );
Lnbr(:, T_BUS)	= e2i( Lnbr(:, T_BUS) );

if size(trsfm,1)>0
   trsfm(:, F_BUS) = e2i( trsfm(:, F_BUS) );
   trsfm(:, T_BUS) = e2i( trsfm(:, T_BUS) );
end

