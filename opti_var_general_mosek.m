
function [loss,gen,bus,Lnbr,trsfm,shtc,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_general_mosek(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)
% 按照MOSEK的SOC约束模型，在P和Q前面需乘以2
% since 2023.02.01

%----------------------------------------------------------------------%
%----------------Pre-Processing-----------------%
% Step 1: Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr);

% Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,~,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr,trsfm,shtc,shtr,sysdt);

nb = size(bus,1); npv=size(pv,1);npq=size(pq,1);
nLnbr = size(Lnbr,1) ;
ng = size(gen,1);
if ~isempty(trsfm)
    ntrsfm = size(trsfm,1);
    frombus = [Lnbr(:,2) ; trsfm(:,2) ];
    tobus = [Lnbr(:,3); trsfm(:,3) ];
else
    frombus = [Lnbr(:,2)  ];
    tobus = [Lnbr(:,3) ];
    ntrsfm = 0;
end

% Step 2:定义优化变量x=[2*P_mn 2*P_nm 2*Q_mn 2*Q_nm I^2_mn I^2_nm U^2 Q_G Qshtc].
% P,Q,I^2 is arranged in the sequence of P_index; U is arranged in the bus index
P_index = [frombus tobus];
A_P = zeros(nb,6*(nLnbr+ ntrsfm)+nb);
delta_P_s = zeros(nb,1);
for i = 2 : nb
    
    % nodes with power flow injections are ones
    ids =find(P_index(:,2)==i);
    A_P(i,(nLnbr+ ntrsfm)+ids) = -1;
    
    %  nodes with power flow ejections are negative ones
    idsm =find(P_index(:,1)==i);
    A_P(i,idsm) = -1;
    
    % nodes with power flow injections are positive ones
    if bus(i,2)==2
        id_g =find(gen(:,1)==i);
        delta_P_s(i,1) = -bus(i,3) - gen(id_g,2);
    else
        delta_P_s(i,1) = -bus(i,3) ;
    end
    %--------end-----------%
    
end
A_P_reduced = 1/2.*A_P;

A_Q = zeros(nb,6*(nLnbr+ ntrsfm)+nb);
delta_Q_s = zeros(nb,1);
for i =  1 : nb
    % ----nodal power balance for node n--------%
    % nodes with power flow ejections are negative ones
    % nodes with power flow injections are ones
    ids =find(P_index(:,2)==i);
    A_Q(i,3*(nLnbr+ ntrsfm)+ids) = -1;
    
    %  nodes with power flow ejections are negative ones
    idsm =find(P_index(:,1)==i);
    A_Q(i,2*(nLnbr+ ntrsfm)+idsm) = -1;
    
    % nodes with power flow injections are positive ones
    delta_Q_s(i,1) = -bus(i,4);
    %--------end-----------%
end

% construct the coincidence matrix of Q_G for all pv buses
npv = size(pv,1);
A_QG = zeros(nb,ng);
for i = 1 : ng
        A_QG( gen(i,1) ,i) = 1;
end

% construct the coincidence shtc matrix
nshtc=size(shtc,1);
A_SHTC = zeros(nb,nshtc);
for i = 1 :nshtc
    A_SHTC( shtc(i,1) ,i) = 1;
end
A_Q_reduced = [1/2.*A_Q A_QG A_SHTC];

% express alpha and constant real coefficient
j = sqrt(-1);
if ~isempty(trsfm)
    alpha_mn = ones(nLnbr+ ntrsfm,1) + ([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]).*[j.*Lnbr(:,6)./2;((1-trsfm(:,6))./(trsfm(:,6).^2)).*1./(trsfm(:,4)+j.*trsfm(:,5))];
    alpha_nm = ones(nLnbr+ ntrsfm,1) + ([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]).*[j.*Lnbr(:,6)./2; (trsfm(:,6)-1)./trsfm(:,6).*1./(trsfm(:,4)+j.*trsfm(:,5))];
    
    const_coef_1 = real(conj(alpha_mn).*[Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))] + alpha_mn.*conj([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]));
    const_coef_2 = imag(conj(alpha_mn).*[Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))] - alpha_mn.*conj([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]));
    
    const_coef_3 = real(conj(alpha_nm).*[Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))] + alpha_nm.*conj([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]));
    const_coef_4 = imag(conj(alpha_nm).*[Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))] - alpha_nm.*conj([Lnbr(:,4)+j.*Lnbr(:,5);trsfm(:,6).*(trsfm(:,4)+j.*trsfm(:,5))]));
    
    z_mn = [(Lnbr(:,5).^2 + Lnbr(:,4).^2);(trsfm(:,6).*trsfm(:,4)).^2 + (trsfm(:,6).*trsfm(:,5)).^2];
    r_mn = [Lnbr(:,4)./2 ; trsfm(:,6).*trsfm(:,4)./2];
    x_mn = [Lnbr(:,5)./2 ; trsfm(:,6).*trsfm(:,5)./2];
else
    alpha_mn = ones(nLnbr+ ntrsfm,1) + ([Lnbr(:,4)+j.*Lnbr(:,5)]).*j.*Lnbr(:,6)./2;
    alpha_nm = ones(nLnbr+ ntrsfm,1) + ([Lnbr(:,4)+j.*Lnbr(:,5)]).*j.*Lnbr(:,6)./2;
    
    const_coef_1 = real(conj(alpha_mn).*[Lnbr(:,4)+j.*Lnbr(:,5)] + alpha_mn.*conj([Lnbr(:,4)+j.*Lnbr(:,5)]));
    const_coef_2 = imag(conj(alpha_mn).*[Lnbr(:,4)+j.*Lnbr(:,5)] - alpha_mn.*conj([Lnbr(:,4)+j.*Lnbr(:,5)]));
    
    const_coef_3 = real(conj(alpha_nm).*[Lnbr(:,4)+j.*Lnbr(:,5)] + alpha_nm.*conj([Lnbr(:,4)+j.*Lnbr(:,5)]));
    const_coef_4 = imag(conj(alpha_nm).*[Lnbr(:,4)+j.*Lnbr(:,5)] - alpha_nm.*conj([Lnbr(:,4)+j.*Lnbr(:,5)]));
    
    z_mn = [(Lnbr(:,5).^2 + Lnbr(:,4).^2)];
    r_mn = [Lnbr(:,4)./2 ;];
    x_mn = [Lnbr(:,5)./2 ;];
end


U_P_var_1 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : nLnbr + ntrsfm
    
    % ----voltage drop for line (m,n)--------%
    U_P_var_1(i,6*(nLnbr+ ntrsfm)+P_index(i,2))  = 1;
    U_P_var_1(i,6*(nLnbr+ ntrsfm)+P_index(i,1))  = -abs(alpha_mn(i,1))^2;
    
    U_P_var_1(i,i) = const_coef_1(i,1)/2;
    U_P_var_1(i,i+2*(nLnbr+ ntrsfm)) = const_coef_2(i,1)/2;
    
    U_P_var_1(i,i+2*(nLnbr+ ntrsfm)+2*(nLnbr+ ntrsfm)) =  -z_mn(i,1);
    % ----end--------%
    
end

U_P_var_2 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    
    % ----voltage drop for line (n,m)--------%
    U_P_var_2(i,6*(nLnbr+ ntrsfm)+P_index(i,1))  = 1;
    U_P_var_2(i,6*(nLnbr+ ntrsfm)+P_index(i,2))  = -abs(alpha_nm(i,1))^2;
    
    U_P_var_2(i,(nLnbr+ ntrsfm)+i) = const_coef_3(i,1)/2;
    U_P_var_2(i,i+3*(nLnbr+ ntrsfm)) = const_coef_4(i,1)/2;
    
    U_P_var_2(i,i+(nLnbr+ ntrsfm)+4*(nLnbr+ ntrsfm)) = -z_mn(i,1);
    % ----end--------%
    
end

br_eq_1 = zeros(nLnbr,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    
    br_eq_1(i,6*(nLnbr+ ntrsfm)+P_index(i,1)) = real(alpha_mn(i,1));
    br_eq_1(i,6*(nLnbr+ ntrsfm)+P_index(i,2)) = -real(alpha_nm(i,1));
    
    br_eq_1(i,i) = -r_mn(i,1);
    br_eq_1(i,2*(nLnbr+ ntrsfm)+i) = -x_mn(i,1);
    
    br_eq_1(i,(nLnbr+ ntrsfm)+i) = r_mn(i,1);
    br_eq_1(i,3*(nLnbr+ ntrsfm)+i) = x_mn(i,1);
    
end

br_eq_2 = zeros(nLnbr,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    
    br_eq_2(i,6*(nLnbr+ ntrsfm)+P_index(i,1)) = -imag(alpha_nm(i,1));
    br_eq_2(i,6*(nLnbr+ ntrsfm)+P_index(i,2)) = -imag(alpha_nm(i,1));
    
    br_eq_2(i,i) = x_mn(i,1);
    br_eq_2(i,2*(nLnbr+ ntrsfm)+i) = -r_mn(i,1);
    
    br_eq_2(i,(nLnbr+ ntrsfm)+i) =  x_mn(i,1);
    br_eq_2(i,3*(nLnbr+ ntrsfm)+i) = -r_mn(i,1);
    
end
% ------------- END ------------------ %

% ------------形成Aeq和Beq矩阵-------- %
% Add auxiliary variables for SOC Constraints
% auxiliary variable W equality constraints
% w1
A_auxiliary_1 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    A_auxiliary_1(i,4*(nLnbr+ ntrsfm) + i) = 1;
    A_auxiliary_1(i,6*(nLnbr+ ntrsfm)+P_index(i,1)) =1;
    A_auxiliary_1(i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+i) = -1;
end
% m1
A_auxiliary_2 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    A_auxiliary_2(i,4*(nLnbr+ ntrsfm) + i) = 1;
    A_auxiliary_2(i,6*(nLnbr+ ntrsfm)+P_index(i,1)) = -1;
    A_auxiliary_2(i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+(nLnbr+ ntrsfm)+i) = -1;
end

% w2
A_auxiliary_3 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    A_auxiliary_3(i,4*(nLnbr+ ntrsfm) +(nLnbr+ ntrsfm)+ i) = 1;
    A_auxiliary_3(i,6*(nLnbr+ ntrsfm)+P_index(i,2)) =1;
    A_auxiliary_3(i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+2*(nLnbr+ ntrsfm)+i) = -1;
end
% m2
A_auxiliary_4 = zeros(nLnbr+ ntrsfm,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm));
for i = 1 : (nLnbr+ ntrsfm)
    A_auxiliary_4(i,4*(nLnbr+ ntrsfm) +(nLnbr+ ntrsfm)+ i) = 1;
    A_auxiliary_4(i,6*(nLnbr+ ntrsfm)+P_index(i,2)) = -1;
    A_auxiliary_4(i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+3*(nLnbr+ ntrsfm)+i) = -1;
end

% Step 3: Form the constant matrix for Aeq*x=beq
Aeq_modify = [ A_P_reduced   zeros(nb,ng+nshtc+4*(nLnbr+ ntrsfm))  ;
    A_Q_reduced  zeros(nb,4*(nLnbr+ ntrsfm))      ;
    U_P_var_1 ;
    U_P_var_2 ;
    br_eq_1;
    br_eq_2;
    A_auxiliary_1  ;
    A_auxiliary_2 ;
    A_auxiliary_3  ;
    A_auxiliary_4 ;
    ];

Beq_modify = [delta_P_s;delta_Q_s;zeros(2*(nLnbr+ ntrsfm),1);zeros(2*(nLnbr+ ntrsfm),1);zeros(4*(nLnbr+ ntrsfm),1);];
% ------------- END ------------------ %

% --------------调用MOSEK求解工具箱-------------%
% Specify the cones.
prob.cones   = cell(2*(nLnbr+ ntrsfm),1);
for i = 1 : nLnbr+ ntrsfm
    prob.cones{i}.type = 'MSK_CT_QUAD';%MSK_CT_RQUAD for Qudratic Cone
    prob.cones{i}.sub  = [6*(nLnbr+ ntrsfm)+nb+ng+nshtc+i,  i,2*(nLnbr+ ntrsfm)+i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+(nLnbr+ ntrsfm)+i];
end

for i = 1 : nLnbr+ ntrsfm
    prob.cones{nLnbr+ ntrsfm+i}.type = 'MSK_CT_QUAD';%MSK_CT_RQUAD for Qudratic Cone
    prob.cones{nLnbr+ntrsfm+i}.sub  = [6*(nLnbr+ ntrsfm)+nb+ng+nshtc+2*(nLnbr+ ntrsfm)+i,  (nLnbr+ ntrsfm)+i,3*(nLnbr+ ntrsfm)+i,6*(nLnbr+ ntrsfm)+nb+ng+nshtc+2*(nLnbr+ ntrsfm)+(nLnbr+ ntrsfm)+i];
end

% Bounds
% 包含U^2电压安全约束、 QG(PV节点的QG)上下限约束、 Qshtc上下限约束
lb = [-10.*ones(1,4*(nLnbr+ ntrsfm)) zeros(1,2*(nLnbr+ ntrsfm)) (bus(:,9).^2)' gen(1:end,5)' shtc(:,4)'  -10.*ones(1,4*(nLnbr+ ntrsfm))];
ub =[10.*ones(1,4*(nLnbr+ ntrsfm))  10.*ones(1,2*(nLnbr+ ntrsfm))  (bus(:,10).^2)' gen(1:end,4)'  shtc(:,3)' 10.*ones(1,4*(nLnbr+ ntrsfm))  ];


% Linear Constraints
% QG(SW节点的QG)上下限约束
% note that the first line contains the PCC bus
% A_QGSW_inequal =  [zeros(1,2*(nLnbr+ ntrsfm)) 1 zeros(1,2*(nLnbr+ ntrsfm)-1) zeros(1,2*(nLnbr+ ntrsfm)) zeros(1,nb) zeros(1,npv+nshtc) zeros(1,4*(nLnbr+ ntrsfm))];
% b_QGSW_inequal_upper = gen(1,4);
% b_QGSW_inequal_lower = gen(1,5);

A =  Aeq_modify;
rl =  Beq_modify;
ru =  Beq_modify ;


% Specify the non-conic part of the problem.
prob.c = [ones(1,2*(nLnbr+ ntrsfm)) zeros(1,4*(nLnbr+ ntrsfm)+nb+ng+nshtc+4*(nLnbr+ ntrsfm)) ];
prob.a = sparse(A);
prob.blc = rl;
prob.buc = ru;
prob.blx = lb;
prob.bux = ub;

% Optimize the problem.
[r,res]=mosekopt('minimize info',prob);
% Display the primal solution.
% res.info.MSK_DINF_INTPNT_PRIMAL_OBJ/2  % 目标函数值
%-------------------------END ------------------------%

X_Opt = res.sol.itr.xx;

% 输出无功优化结果
busV = sqrt(X_Opt(6*(nLnbr+ ntrsfm)+1:6*(nLnbr+ ntrsfm)+nb));  % 最优电压幅值
bus(:,7) = busV;
bus(:,8) = 0;

ids = find(P_index(:,1)==1);
idm = find(P_index(:,2)==1);
Pns = (sum(X_Opt(ids ,1),1) + sum(X_Opt(nLnbr+ ntrsfm+idm ,1),1))/2;% 最优注入有功

gen(1,2) = Pns;
gen(1,3) = (sum(X_Opt(2*(nLnbr+ ntrsfm) + ids ,1),1) + sum(X_Opt(2*(nLnbr+ ntrsfm) + nLnbr+ ntrsfm+idm ,1),1))/2;

gen(1,6) = sqrt(X_Opt(6*(nLnbr+ ntrsfm)+1,1));
for i =1 : ng
    gen(i,3) = X_Opt(6*(nLnbr+ ntrsfm)+nb+i,1);
    gen(i,6) = busV(gen(i,1),1);
end

shtc(:,2) = X_Opt(6*(nLnbr+ ntrsfm)+nb+ng+1:6*(nLnbr+ ntrsfm)+nb+ng+nshtc,1);   % 最优注入无功补偿
loss = res.info.MSK_DINF_INTPNT_PRIMAL_OBJ/2;  % 最小网损

Lnbr(:,7) = [X_Opt(1:nLnbr)./2]';
Lnbr(:,8) = [X_Opt(2*(nLnbr+ ntrsfm)+1:2*(nLnbr+ ntrsfm)+nLnbr)./2]';
Lnbr(:,9) = X_Opt((nLnbr+ ntrsfm)+1:(nLnbr+ ntrsfm)+nLnbr)./2;
Lnbr(:,10) =  X_Opt(3*(nLnbr+ ntrsfm)+1:3*(nLnbr+ ntrsfm)+nLnbr)./2;
Lnbr(:,11) = bus(Lnbr(:,2),7);
Lnbr(:,12) = bus(Lnbr(:,3),7);

if ~isempty(trsfm)
    trsfm(:,11) = [X_Opt(1+nLnbr:nLnbr+ ntrsfm)./2]';
    trsfm(:,12) = [X_Opt(2*(nLnbr+ ntrsfm)+nLnbr+1:3*(nLnbr+ ntrsfm))./2]';
    trsfm(:,13) = X_Opt((nLnbr+ ntrsfm)+nLnbr+1:2*(nLnbr+ ntrsfm))./2;
    trsfm(:,14) = X_Opt(3*(nLnbr+ ntrsfm)+nLnbr+1:4*(nLnbr+ ntrsfm))./2;
    trsfm(:,15) = bus(trsfm(:,2),7);
    trsfm(:,16) = bus(trsfm(:,3),7);
end


% Step 4: Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if res.rcode == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

% Nonlinear Constraints（SOC约束）
nlcon_1 = @(x) [ x( (nLnbr+ ntrsfm)*6+nb+ng+nshtc+1:(nLnbr+ ntrsfm)*6+nb+ng+nshtc+(nLnbr+ ntrsfm) ).^2  ...
    -  x( 1: (nLnbr+ ntrsfm)).^2 - ( x(2*(nLnbr+ ntrsfm) + 1:(nLnbr+ ntrsfm)*2 + (nLnbr+ ntrsfm)) ).^2 ...
    - (x( (nLnbr+ ntrsfm)*6+nb+ng+nshtc+(nLnbr+ ntrsfm)+1:(nLnbr+ ntrsfm)*6+nb+ng+nshtc+2*(nLnbr+ ntrsfm) ).^2) ]';

nlcon_2 = @(x) [ x((nLnbr+ ntrsfm)*6+nb+ng+nshtc+2*(nLnbr+ ntrsfm)+1:(nLnbr+ ntrsfm)*6+nb+ng+nshtc+3*(nLnbr+ ntrsfm) ).^2  ...
    - x( (nLnbr+ ntrsfm)+1: 2*(nLnbr+ ntrsfm)).^2 - x(3*(nLnbr+ ntrsfm) + 1:(nLnbr+ ntrsfm)*3 + (nLnbr+ ntrsfm)).^2 ...
    - x((nLnbr+ ntrsfm)*6+nb+ng+nshtc+3*(nLnbr+ ntrsfm)+1:(nLnbr+ ntrsfm)*6+nb+ng+nshtc+4*(nLnbr+ ntrsfm)).^2 ]';

SOC_value = max([nlcon_1(X_Opt),nlcon_2(X_Opt)]);


I2 = [X_Opt((nLnbr+ ntrsfm)*4 +1 :(nLnbr+ ntrsfm)*5,1),X_Opt((nLnbr+ ntrsfm)*5 +1 :(nLnbr+ ntrsfm)*6,1)];


sysdt(PFITER) = res.info.MSK_IINF_INTPNT_ITER; % the number of interior-point iterations

% Step 6: Bus Solution and Line Flows
[gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = int2ext(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss);

elapsed_time = res.info.MSK_DINF_OPTIMIZER_TIME;  %  total optimization time
