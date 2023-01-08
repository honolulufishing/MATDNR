
function [loss,gen,bus,Lnbr,trsfm,shtc,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_mosek(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)
% 按照MOSEK的SOC约束模型，在P和Q前面需乘以2
% since 2020.03.22

%----------------Pre-Processing-----------------%
% Step 1: Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr);

% Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr,trsfm,shtc,shtr,sysdt);
 
frombus = Lnbr(:,2);
tobus=Lnbr(:,3);
nb = size(bus,1); npv=size(pv,1);npq=size(pq,1);
nLnbr=nb-1;

% Step 2:定义优化变量[sqrt(2)*P sqrt(2)*Q I^2 U^2 QGforPV Qshtc].其中，P,Q,I^2,U^2 is in the dimensions of (nb-1)*1,(nb-1)*1,(nb-1)*1,and (nb)*1.
% P,Q,I^2 is arranged in the sequence of P_index; U is arranged in the bus
% index However, note that 点2不能是PV节点.
idgen=zeros(npv,1);
for j = 1: npv
idgen(j,1) =find(gen(:,1)==pv(j,1));
end

P_index = [frombus tobus];
A_P = zeros(nb,nLnbr+2*nLnbr);
for i = 2 : nb
    % nodes with power flow injections are ones
    ids =find(P_index(:,2)==i);
    A_P(i,ids) = 1;
    A_P(i,2*nLnbr+ids) = -2*Lnbr(ids,4);
    %  nodes with power flow ejections are negative ones
    idsm =find(P_index(:,1)==i);
    A_P(i,idsm) = -1;
end
A_P_reduced = A_P(2:nb,1:3*nLnbr);

A_Q = zeros(nb,nLnbr+2*nLnbr);
for i = 2 : nb
    % nodes with power flow injections are ones
    ids =find(P_index(:,2)==i);
    A_Q(i,nLnbr+ids) = 1;
    A_Q(i,2*nLnbr+ids) = -2*Lnbr(ids,5);
    %  nodes with power flow ejections are negative ones
    idsm =find(P_index(:,1)==i);
    A_Q(i,nLnbr+idsm) = -1;
end
A_Q_reduced = A_Q(2:nb,1:3*nLnbr);


U_P_var = zeros(nb-1,nb-1+nb-1+nb-1+nb);
for i = 1 :nb-1

       U_P_var(i,nb-1+nb-1+nb-1+P_index(i,1))  = 1;  %根据节点序号来设置电压变量顺序
       U_P_var(i,nb-1+nb-1+nb-1+ P_index(i,2))  = -1;
        
        U_P_var(i,i) = -Lnbr(i,4);
        U_P_var(i,i+nb-1) = -Lnbr(i,5);
        U_P_var(i,i+(nb-1)*2) = Lnbr(i,5)^2 + Lnbr(i,4)^2;

end
% ------------- END ------------------ %

%----------------采用MOSEK求解无功优化-----------------%
% x=[2P 2Q I^2 U^2 QG Qshtc]:
% U^2电压安全约束
% QG上下限约束
% Qshtc上下限约束

% 添加QG矩阵（PV节点）的顺序来自PV节点顺序
A_UQG = zeros(nb-1,npv);
for i =1:npv
    idf = find(P_index(:,2)==pv(i));
    if ~isempty(idf)
        A_UQG(idf,i) =1;
    end
end

% construct the coincidence shtc matrix
nshtc=size(shtc,1);
A_SHTC = zeros(nb-1,nshtc);
for i = 1 :nshtc
    A_SHTC(shtc(i,1)-1,i) =1;
end


delta_P_s = -bus(2:end,3) - [zeros(npq,1);gen(idgen,2)];
delta_Q_s = -bus(2:end,4) ;

% Add auxiliary variables for SOC Constraints
% auxiliary variable W equality constraints
A_auxiliary_1 = zeros(nb-1,(nb-1)+nb+(nb-1)*2+npv+nshtc);
for i = 1 :nb-1
    A_auxiliary_1(i,i) = 1;
    A_auxiliary_1(i,nb-1+P_index(i,1)) =1;
     A_auxiliary_1(i,(nb-1)+nb+npv+nshtc+i) = -1;
end
% m
A_auxiliary_2 = zeros(nb-1,(nb-1)+nb+(nb-1)*2+npv+nshtc);
for i = 1 :nb-1
    A_auxiliary_2(i,i) = 1;
    A_auxiliary_2(i,nb-1+P_index(i,1)) =-1;
      A_auxiliary_2(i,(nb-1)+nb+npv+nshtc+nb-1+i) = -1;
end

% ------------形成Aeq和Beq矩阵-------- %
% Step 3: Form the constant matrix for Aeq*x=beq
Aeq_modify = [ (1/2).*A_P_reduced   zeros(nb-1,nb+npv+nshtc) zeros(nb-1,(nb-1)*2) ;
   (1/2).*A_Q_reduced   zeros(nb-1,nb)   A_UQG A_SHTC zeros(nb-1,(nb-1)*2); 
        U_P_var zeros(nb-1,npv+nshtc) zeros(nb-1,(nb-1)*2)
          zeros(nb-1,nb-1+nb-1)  A_auxiliary_1  ;
         zeros(nb-1,nb-1+nb-1) A_auxiliary_2 ];

Beq_modify = [delta_P_s;delta_Q_s;zeros(nb-1,1);zeros(2*(nb-1),1)];
% ------------- END ------------------ %

% --------------调用MOSEK求解工具箱-------------%
% Specify the cones.
prob.cones   = cell(nb-1,1);
for i = 1 : nb-1
    prob.cones{i}.type = 'MSK_CT_QUAD';%MSK_CT_RQUAD for Qudratic Cone
    prob.cones{i}.sub  = [3*(nb-1)+nb+nshtc+npv+i,  i,nb-1+i,3*(nb-1)+nb+nshtc+npv+(nb-1)+i];
end

% Bounds
% 包含U^2电压安全约束、 QG(PV节点的QG)上下限约束、 Qshtc上下限约束
lb = [-100.*ones(1,3*(nb-1)) (bus(:,9).^2)' gen(2:end,5)'  shtc(:,4)' -100.*ones(1,2*(nb-1))  ];
ub =[100.*ones(1,3*(nb-1)) (bus(:,10).^2)' gen(2:end,4)'  shtc(:,3)'   100.*ones(1,2*(nb-1))];


% Linear Constraints
% QG(SW节点的QG)上下限约束
A_QGSW_inequal =  [zeros(1,nb-1) 1 zeros(1,nb-2) zeros(1,nb-1) zeros(1,nb) zeros(1,npv) zeros(1,nshtc) zeros(1,2*(nb-1))];
b_QGSW_inequal_upper = gen(1,4);
b_QGSW_inequal_lower = gen(1,5);

A = [ Aeq_modify
    A_QGSW_inequal ];
rl = [ Beq_modify;b_QGSW_inequal_lower];
ru = [Beq_modify ;b_QGSW_inequal_upper];

% Specify the non-conic part of the problem.
prob.c = [1 zeros(1,(nb-1)*3+nb +nshtc+npv-1+2*(nb-1)) ];
prob.a = sparse(A);
prob.blc = rl;
prob.buc = ru;
prob.blx = lb;
prob.bux = ub;

% Optimize the problem.
[r,res]=mosekopt('minimize info',prob);
% Display the primal solution.
% res.info.MSK_DINF_INTPNT_PRIMAL_OBJ  % 目标函数值
%-------------------------END ------------------------%

X_Opt = res.sol.itr.xx;

% 输出无功优化结果
busV = sqrt(X_Opt(3*(nb-1)+1:3*(nb-1)+nb));  % 最优电压幅值
bus(:,7) = busV;
bus(:,8) = 0;

Pns = X_Opt(1,1)./2;                   % 最优注入有功
gen(pv,3) = X_Opt(3*(nb-1)+nb+1:3*(nb-1)+nb+npv);
gen(1,3) = X_Opt((nb-1)+1)./2;
gen(1,2) = Pns;
gen(1,6) = sqrt(X_Opt(3*(nb-1)+1,1));


shtc(:,2) = X_Opt(3*(nb-1)+nb+1,1);           % 最优注入无功补偿
loss = X_Opt(1,1)./2+sum(bus(:,3),1);  % 最小网损

Lnbr(:,7) = [X_Opt(1:nb-1)./2]';
Lnbr(:,8) = [X_Opt((nb-1)+1:(nb-1)+nb-1)./2]';

Lnbr(:,9) = - (X_Opt(1:nb-1)./2 - [X_Opt((nb-1)*2+1:(nb-1)*2+nb-1)].*Lnbr(:,4));
Lnbr(:,10) = -(X_Opt((nb-1)+1:(nb-1)+nb-1)./2 - [X_Opt((nb-1)*2+1:(nb-1)*2+nb-1)].*Lnbr(:,5));

Lnbr(:,11) = bus(Lnbr(:,2),7);
Lnbr(:,12) = bus(Lnbr(:,3),7);

% Step 4: Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if res.rcode == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

% Nonlinear Constraints（SOC约束）
nlcon = @(x) [ x( (nb-1)*3+nb +npv+nshtc+1:(nb-1)*3+nb +npv+nshtc+(nb-1) ).^2  ...
    - ( x( 1: nb-1) ).^2 - ( x(nb-1 + 1:(nb-1)*1 + nb-1) ).^2 ...
    - (x( (nb-1)*3+nb+npv+nshtc+nb-1+1:(nb-1)*3+nb +npv+nshtc+2*(nb-1) ).^2) ]';     

SOC_value = nlcon(X_Opt);


I2 =X_Opt((nb-1)*2 +1 :(nb-1)*2 + nb-1);


sysdt(PFITER) = res.info.MSK_IINF_INTPNT_ITER; % the number of interior-point iterations

% Step 6: Bus Solution and Line Flows 
[gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = int2ext(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss);

elapsed_time = res.info.MSK_DINF_OPTIMIZER_TIME;  %  total optimization time


