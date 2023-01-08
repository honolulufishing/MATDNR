function [gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = int2ext_NR(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss,beta)
% Converts internal to external bus numbering.

% define names for columns to data matrices
[GEN_BUS,GEN_PG,GEN_QG,QGMAX,QGMIN,GEN_VG,GEN_OPTVG,GEN_OPTQG,GEN_OPTPG] = idx_Gen;
[LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,LN_OPTPF,LN_OPTQF,...
      LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR;
[BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus;
[TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,TRS_PF,TRS_QF,TRS_PT,TRS_QT,...
      TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM;
[VCTR_BUS,VSET,WEIGHT,VCTR_VM,VCTR_OPTVM] = idx_VCTR;
[SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC;
[SHTR_BUS,SHTR_QR,QRMAX,QRMIN,SHTR_SERIES,SHTR_STATE,OPTQR] = idx_SHTR;

bus(:, BUS_I) = i2e( bus(:, BUS_I) );
gen(:, GEN_BUS)	= i2e( gen(:, GEN_BUS) );
if size(shtc,1)>0, shtc(:,SHTC_BUS)	= i2e( shtc(:, SHTC_BUS) );end
if size(shtr,1)>0, shtr(:,SHTR_BUS)	= i2e( shtr(:, SHTR_BUS) );end
if size(vctr,1)>0, vctr(:,VCTR_BUS)	= i2e( vctr(:, VCTR_BUS) );end

Lnbr(:, F_BUS)	= i2e( Lnbr(:, F_BUS) );
Lnbr(:, T_BUS)	= i2e( Lnbr(:, T_BUS) );
if size(trsfm,1)>0
   trsfm(:, F_BUS)	= i2e( trsfm(:, F_BUS) );
   trsfm(:, T_BUS)	= i2e( trsfm(:, T_BUS) );
end

% Output Solution
format short
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

% Output Bus Solution
disp('                                                                            ')
disp('                  Optimal DNR Solution by MATDNR Toolbox')
fprintf('                 The System Total Real Loss = %g MW\n\n', loss)
head =['    Bus  Bus    ------Load------    ---Injected---     Voltage    Angle  '
       '    No.  Type    MW       Mvar       GS       BS         Mag.      Deg.  '
       '                                                                         '];
disp(head)
for i = 1:size(bus,1)
    fprintf(' %5g', bus(i,1)), fprintf(' %5g', bus(i,2)),
    fprintf(' %8.3f', bus(i,3)), fprintf(' %9.3f', bus(i,4)),
    fprintf(' %9.3f', bus(i,5)),  fprintf(' %9.3f', bus(i,6)),
    fprintf(' %9.3f ', bus(i,7)), fprintf(' %8.3f\n',bus(i,8)*180/pi)
end
fprintf('      \n'), fprintf('  Total     ')
   fprintf(' %8.3f', sum(bus(:,3))), fprintf(' %9.3f\n\n', sum(bus(:,4)))


% Output Generator Solution
disp('                                     ')
disp('                   Generator Solution')

head =['                                                     '
       '    Bus     ----Gen----       -----Q-----       Vol  '
       '    No.     MW      MVA       Max     Min       Mag. '
       '                                                     '];
disp(head)
for i = 1:size(gen,1)
    fprintf(' %5g', gen(i,1)), fprintf(' %8.3f', gen(i,2)), 
    fprintf(' %8.3f', gen(i,3)), fprintf(' %9.3f',gen(i,4)),
    fprintf(' %9.3f', gen(i,5)),  fprintf(' %8.3f\n', gen(i,6)),
end
fprintf('      \n'), fprintf('  Total ')
fprintf(' %6.3f', sum(gen(:,2))), fprintf(' %8.3f\n\n', sum(gen(:,3)))
   
% Output Line Solution
disp('                                                  ')
disp('                   Transmission Line Flow and Loss')

head =['                                                                                                                                            '
       '   LnBR.  Bus   Bus    R      X    1/2 B      Power Flow        Power Flow           Loss          VolMag   VolMag    Switch     --beta--   '
       '    NO.  from   to    p.u.   p.u.   p.u.        f to t            t to f             p.u.            f        t       Status     nm   mn    '
       '                                                                                                                                            '];
disp(head)

for i = 1:size(Lnbr,1)
    fprintf(' %5g', Lnbr(i,1)), fprintf(' %5.0f', Lnbr(i,2)),
    fprintf(' %5g', Lnbr(i,3)), fprintf(' %8.3f', Lnbr(i,4)),
    fprintf(' %6.3f', Lnbr(i,5)),fprintf(' %6.3f', Lnbr(i,6)),
    fprintf(' %8.3f+j%7.3f', Lnbr(i,7),Lnbr(i,8)),
    fprintf(' %8.3f+j%7.3f', Lnbr(i,9),Lnbr(i,10)),
    fprintf(' %8.3f+j%7.3f', Lnbr(i,7)+Lnbr(i,9),Lnbr(i,8)+Lnbr(i,10)),
    fprintf(' %8.3f', Lnbr(i,11)), fprintf(' %8.3f', Lnbr(i,12)),
    fprintf(' %8.3f', Lnbr(i,13)),fprintf(' %8.3f', beta(2*i,1)),fprintf(' %8.3f\n', beta(2*i-1,1)),
end
fprintf('      \n'), fprintf('    Total                                ')
fprintf(' %8.3f+j%7.3f', sum(Lnbr(:,7)),sum(Lnbr(:,8))),
fprintf(' %8.3f+j%7.3f', sum(Lnbr(:,9)),sum(Lnbr(:,10))),
fprintf(' %8.3f+j%7.3f\n\n', sum(Lnbr(:,7)+Lnbr(:,9)),sum(Lnbr(:,8)+Lnbr(:,10)))


% ---- Additonal Facility Data ---- %
% ----------- shunt capacity Data ------------ %
%      Bus  shtc     ---Q---  
%       NO. Mvar  Max  Min 


% Output Additonal Facility Solution
disp('                                                  ')
disp('                   Additonal Facility Shunt Capacity Solution')

head =['                                    '
       '   Bus.   shtc     ----Q----        '
       '    NO.   Mvar     Max    Min       '
       '                                    '];
disp(head) 

for i = 1:size(shtc,1)
    fprintf(' %5g', shtc(i,1)), fprintf('   %5.3f', shtc(i,2)),
    fprintf('   %5g', shtc(i,3)), fprintf('   %5.3f\n', shtc(i,4)),
end

% Output Transformer Branch Solution
disp('                                                  ')
disp('                   Transformer Branch Flow and Loss')
 
head =['                                                                                                                                               '
       '   TrsBR.  Bus   Bus     R      X     ----- Tap -----    Tap     Tap      Power Flow      Power Flow           Loss          VolMag    VolMag  '
       '    NO.    from   to    p.u.   p.u.   Ratio  Max  Min   Series  Status      f to t          t to f              p.u.            f        t     '
       '                                                                                                                                               '];
disp(head)

for i = 1:size(trsfm,1)
    fprintf(' %5g', trsfm(i,1)), fprintf(' %5.0f', trsfm(i,2)),
    fprintf(' %5g', trsfm(i,3)), fprintf(' %8.3f', trsfm(i,4)),
    fprintf(' %6.3f', trsfm(i,5)),fprintf(' %5.1f',trsfm(i,6)),
    fprintf(' %5.1f', trsfm(i,7)),fprintf(' %5.1f',trsfm(i,8)),
    fprintf(' %6.0f',trsfm(i,9)),fprintf(' %8.0f',trsfm(i,10)),
    fprintf(' %8.3f+j%7.3f', trsfm(i,11),trsfm(i,12)),
    fprintf(' %8.3f+j%7.3f', trsfm(i,13),trsfm(i,14)),
    fprintf(' %8.3f+j%7.3f',trsfm(i,11)+trsfm(i,13),trsfm(i,12)+trsfm(i,14)),
    fprintf(' %8.3f', trsfm(i,15)), fprintf(' %8.3f\n', trsfm(i,16)),
end
fprintf('      \n'), fprintf('    Total                                ')
fprintf(' %8.3f+j%7.3f', sum(trsfm(:,11)),sum(trsfm(:,12))),
fprintf(' %8.3f+j%7.3f', sum(trsfm(:,13)),sum(trsfm(:,14))),
fprintf(' %8.3f+j%7.3f\n\n', sum(trsfm(:,11)+trsfm(:,13)),sum(trsfm(:,12)+trsfm(:,14)))








