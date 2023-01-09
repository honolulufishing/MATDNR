
function [loss,gen,bus,Lnbr,trsfm,shtc,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_mosek_new(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)
% 按照MOSEK的SOC约束模型，在P和Q前面需乘以2
% since 2020.03.22

%----------------Pre-Processing-----------------%
% Step 1: Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus_s,Lnbr_s,trsfm_s,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr);

% Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus_s,Lnbr_s,trsfm_s,shtc,shtr,sysdt);
% --------------------- END -------------------- %

% ----------- 构造变压器支路的虚拟节点 ---------- %
nb = size(bus,1);

    ntb = size(find(trsfm(:,10)==1),1); % 新增虚拟节点个数
    ntrsfm = ntb;
    
    temp = find(trsfm(:,10) == 1 );
    tbnode = [trsfm(temp,2),trsfm(temp,3) [nb+1:nb+ntb]' trsfm(temp,3) ];
    trsfm(temp,2) = [nb+1:nb+ntb]'; % 变压器节点变化转为虚拟节点编号
    
    % 引入虚拟节点
    bus(nb+1:nb+ntb,1) = tbnode(:,3);
    bus(nb+1:nb+ntb,8) = bus(tbnode(:,3),8);
    bus(nb+1:nb+ntb,2) = 1;  % 虚拟节点为PQ节点类型
    
    bus(find(bus(:,9)==0),9)= 0.8;     % 虚拟节点的电压下限
    bus(find(bus(:,10)==0),10)= 1.2 ;  % 虚拟节点的电压上限
    
    Lnbr_temp = [trsfm(:,1:5) 0 ];
    Lnbr= [Lnbr_temp ;Lnbr];
    Lnbr(:,1) = [1:nb-1]';  %nb为真实节点个数（不包含虚拟节点）


% --------------------- END -------------------- %

%----------------Pre-Processing continued-----------------%
frombus = Lnbr(:,2);
tobus=Lnbr(:,3);
npv=size(pv,1);npq=size(pq,1);
NB_modify = size(bus,1);  % number of all buses including added virtual buses

% Step 2:定义优化变量[2*P 2*Q I^2 U^2 QGforPV Qshtc 1/2tap^2].其中，P,Q,I^2,U^2 is in the dimensions of (nb-1)*1,(nb-1)*1,(nb-1)*1,and (nb)*1.
% P,Q,I^2 is arranged in the sequence of P_index; U is arranged in the bus
% index However, note that 点2不能是PV节点.
P_index = [frombus tobus];

idgen=zeros(npv,1);
for j = 1: npv
idgen(j,1) =find(gen(:,1)==pv(j,1));
end

% 节点2至Nb的潮流方程by[pq;pv]
A_P = zeros(nb-1,nb-1);
for i = 1 : nb-1
    
    %A_P流入节点的为1的赋值
    ids =find(P_index(:,2)==i+1);
    A_P(i,ids) = 1;
    
    %A_P流出节点的为-1的赋值
    idsm =find(P_index(:,1)==i+1);
    A_P(i,idsm) =- 1;
    
end
A_Q = A_P;

I_P =zeros(nb-1,nb-1);
I_Q =zeros(nb-1,nb-1);
for i = 1 :nb-1
    %I_P流入节点的为-1的赋值
    ids =find(P_index(:,2)==i+1);
    
    tf = find(frombus==P_index( ids ,1));
    tt = find(tobus(tf)==P_index(ids ,2));
    if ~isempty(tt)
        I_P(i,ids) = -Lnbr(tf(tt),4);
        I_Q(i,ids) = -Lnbr(tf(tt),5);
    end
    
end

% develop U_P矩阵
U_P_var = zeros(nb-1,nb-1+nb-1+nb-1+NB_modify);
for i = 1 :nb-1
    
    %U_P流入节点的为1的赋值
    ids =find(P_index(:,2)==i+1);
    tf = find(frombus==P_index( ids ,1));
    tt = find(tobus(tf)==P_index(ids ,2));
  
    if ~isempty(tt)
       U_P_var(i,nb-1+nb-1+nb-1+P_index(ids,1))  = 1;  %根据节点序号来设置电压变量顺序
       U_P_var(i,nb-1+nb-1+nb-1+ P_index(ids,2))  = -1;
        
        U_P_var(i,ids) = -Lnbr(tf(tt),4);
        U_P_var(i,ids+nb-1) = -Lnbr(tf(tt),5);
        U_P_var(i,ids+(nb-1)*2) = Lnbr(tf(tt),5)^2 + Lnbr(tf(tt),4)^2;
    end
 
end

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

% 添加Shtc矩阵的顺序来自shtc节点顺序
nshtc=size(shtc,1);
A_SHTC = zeros(nb-1,nshtc);
for i =1:nshtc
    idf = find(P_index(:,2)==shtc(i,1));
    if ~isempty(idf)
        A_SHTC (idf,i) =1;
    end
end


delta_P_s = -bus(2:nb,3) - [zeros(npq,1);gen(idgen,2)];
delta_Q_s = -bus(2:nb,4) ;
% ------------- END ------------------ %

%----------------采用MOSEK求解无功优化-----------------%
% Add auxiliary variables for SOC Constraints
% auxiliary variable W equality constraints
A_auxiliary_1 = zeros(nb-1,(nb-1)+NB_modify+(nb-1)*2+npv+nshtc+ntrsfm);
for i = 1 :nb-1
    A_auxiliary_1(i,i) = 1;
    A_auxiliary_1(i,nb-1+P_index(i,1)) =1;
     A_auxiliary_1(i,(nb-1)+NB_modify+npv+nshtc+ntrsfm+i) = -1;
end
% m
A_auxiliary_2 = zeros(nb-1,(nb-1)+NB_modify+(nb-1)*2+npv+nshtc+ntrsfm);
for i = 1 :nb-1
    A_auxiliary_2(i,i) = 1;
    A_auxiliary_2(i,nb-1+P_index(i,1)) =-1;
      A_auxiliary_2(i,(nb-1)+NB_modify+npv+nshtc+ntrsfm+nb-1+i) = -1;
end

% ------------形成Aeq和Beq矩阵-------- %
% Step 3: Form the constant matrix for Aeq*x=beq
Aeq_modify = [ (1/2).*A_P zeros(nb-1,nb-1) I_P zeros(nb-1,NB_modify+npv+nshtc+ntrsfm) zeros(nb-1,(nb-1)*2) ;
       zeros(nb-1,nb-1) (1/2).*A_Q I_Q zeros(nb-1,NB_modify) A_UQG A_SHTC zeros(nb-1,ntrsfm) zeros(nb-1,(nb-1)*2) ; 
        U_P_var zeros(nb-1,npv+nshtc+ntrsfm) zeros(nb-1,(nb-1)*2) 
        zeros(nb-1,nb-1+nb-1)  A_auxiliary_1  ;
         zeros(nb-1,nb-1+nb-1) A_auxiliary_2 ];

Beq_modify = [delta_P_s;delta_Q_s;zeros(nb-1,1);zeros(2*(nb-1),1)];
% ------------- END ------------------ %

% --------------调用MOSEK求解工具箱-------------%
% Step 4: Solve with BARON toolbox
[r, res] = mosekopt('symbcon');

% Specify the cones.
prob.cones   = cell(nb-1+ntrsfm,1);
for i = 1 : nb-1+ntrsfm
    prob.cones{i}.type = 'MSK_CT_QUAD';% Qudratic Cone
    prob.cones{i}.sub  = [3*(nb-1)+NB_modify+nshtc+npv+ntrsfm+i,  i,nb-1+i,3*(nb-1)+NB_modify+nshtc+npv+ntrsfm+(nb-1)+i];
    if i>nb-1
        prob.cones{i}.type = 'MSK_CT_RQUAD'; % Rotated Qudratic Cone
        prob.cones{i}.sub  = [3*(nb-1)+NB_modify+nshtc+npv+i-nb+1,  3*(nb-1)+trsfm(i-nb+1,2),3*(nb-1)+tbnode(i-nb+1,1)];
    end
    
end


% Bounds
% 包含U^2电压安全约束、 QG(PV节点的QG)上下限约束、 Qshtc上下限约束
lb = [-100.*2.*ones(1,3*(nb-1)) 0.9.^2*ones(1,NB_modify) gen(2:end,5)'  shtc(:,4)'  0.5.*trsfm(:,8)' -100.*ones(1,2*(nb-1))  ];
ub =[100.*ones(1,3*(nb-1)) 1.1.^2*ones(1,NB_modify) gen(2:end,4)'  shtc(:,3)'    0.5.*trsfm(:,7)' 100.*ones(1,2*(nb-1))];

% Linear Constraints
% QG(SW节点的QG)上下限约束
A_QGSW_inequal =  [zeros(1,nb-1) 1 zeros(1,nb-2) zeros(1,nb-1) zeros(1,NB_modify) zeros(1,npv) zeros(1,nshtc+ntrsfm) zeros(1,2*(nb-1))];
b_QGSW_inequal_upper = gen(1,4);
b_QGSW_inequal_lower = gen(1,5);

A = [ Aeq_modify
    A_QGSW_inequal ];
rl = [ Beq_modify;b_QGSW_inequal_lower];
ru = [Beq_modify ;b_QGSW_inequal_upper];

% Specify the non-conic part of the problem.
prob.c = [1 zeros(1,(nb-1)*3+NB_modify +nshtc+npv+ntrsfm-1+2*(nb-1)) ];
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

x = res.sol.itr.xx;

% 输出无功优化结果
busV = sqrt(x(3*(nb-1)+1:3*(nb-1)+NB_modify));  % 最优电压幅值
busV(tbnode(:,3),:)=[];  %消去新增节点
bus_s(1:nb,7) = busV;
bus_s(1:nb,8) = 0;

trsfm_s(:,6) = sqrt(2.*x(3*(nb-1)+NB_modify+npv+nshtc+1:3*(nb-1)+NB_modify+npv+nshtc+ntrsfm));


Pns = x(1,1)./2;                   % 最优注入有功
gen(pv,3) = x(3*(nb-1)+nb+1:3*(nb-1)+NB_modify+npv);
gen(1,3) = x((nb-1)+1)./2;
gen(1,2) = Pns;

shtc(:,2) = x(3*(nb-1)+NB_modify+1,1);           % 最优注入无功补偿
loss = x(1,1)./2+sum(bus(:,3),1);  % 最小网损

Lnbr(1:ntrsfm,:)=[]; %复原在Lnbr原始节点编号

Lnbr(:,7) = [x(ntrsfm+1:nb-1)./2]';          
Lnbr(:,8) = [x((nb-1)+ntrsfm+1:(nb-1)+nb-1)./2]';

Lnbr(:,9) = - (x(ntrsfm+1:nb-1)./2 - [x((nb-1)*2+ntrsfm+1:(nb-1)*2+nb-1)].*Lnbr(:,4));
Lnbr(:,10) = -(x((nb-1)+ntrsfm+1:(nb-1)+nb-1)./2 - [x((nb-1)*2+ntrsfm+1:(nb-1)*2+nb-1)].*Lnbr(:,5));

Lnbr(:,11) = bus(Lnbr(:,2),7);  
Lnbr(:,12) = bus(Lnbr(:,3),7);

trsfm_s(:,11) = [x(1:ntrsfm)./2]';          
trsfm_s(:,12) = [x((nb-1)+1:(nb-1)+ntrsfm)./2]';

trsfm_s(:,13) = - (x(1:ntrsfm)./2 - [x((nb-1)*2+1:(nb-1)*2+ntrsfm)].*trsfm_s(:,4));
trsfm_s(:,14) = -(x((nb-1)+1:(nb-1)+ntrsfm)./2 - [x((nb-1)*2+1:(nb-1)*2+ntrsfm)].*trsfm_s(:,5));

trsfm_s(:,15) = bus(trsfm_s(:,2),7);  
trsfm_s(:,16) = bus(trsfm_s(:,3),7);

I2 =x((nb-1)*2 +1 :(nb-1)*2 + nb-1);

% Step 4: Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if res.rcode == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

% Nonlinear Constraints（SOC约束）
nlcon = @(x) [ x( (nb-1)*3+NB_modify +npv+nshtc+ntrsfm+1:(nb-1)*3+NB_modify +npv+nshtc+ntrsfm+(nb-1) ).^2  ...
    - ( x( 1: nb-1) ).^2 - ( x(nb-1 + 1:(nb-1)*1 + nb-1) ).^2 ...
    - (x( (nb-1)*3+NB_modify+npv+nshtc+ntrsfm+nb-1+1:(nb-1)*3+NB_modify +npv+nshtc+ntrsfm+2*(nb-1) ).^2);
    x( (nb-1)*3+NB_modify +npv+nshtc+1:(nb-1)*3+NB_modify +npv+nshtc+ntrsfm ).*x( (nb-1)*3+trsfm(:,2)) ...
    - x( (nb-1)*3+tbnode(:,1))
    ]';   

SOC_value = nlcon(x);

I2 =x((nb-1)*2 +1 :(nb-1)*2 + nb-1);

sysdt(PFITER) = res.info.MSK_IINF_INTPNT_ITER; % the number of interior-point iterations

% Step 6: Bus Solution and Line Flows 
[gen,bus_s,Lnbr,trsfm_s,shtc,shtr,vctr] = int2ext_trs(i2e, gen,bus_s,Lnbr,trsfm_s,shtc,shtr,vctr,sysdt,loss);

elapsed_time = res.info.MSK_DINF_OPTIMIZER_TIME;  %  total optimization time


