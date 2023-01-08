
function [loss,gen,bus,Lnbr,trsfm,shtc,shtr,vctr,elapsed_time,SOC_value,I2] = opti_varSDP_sedumi(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)
% 按照SeduMi的SDP with complex values约束模型，在P和Q前面需乘以2
% since 2020.03.31

%Lnbr= Lnbr_new;
clear prob

%----------------Pre-Processing-----------------%
% Step 1: Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr,trsfm,shtc,shtr,vctr);

% Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr,trsfm,shtc,shtr,sysdt);
 
frombus = Lnbr(:,2);
tobus=Lnbr(:,3);
nb = size(bus,1); npv=size(pv,1);npq=size(pq,1);
nLnbr= nb-1;

% Step 2:定义优化变量[2*P 2*Q I^2 U^2 QGforPV Qshtc].其中，P,Q,I^2,U^2 is in the dimensions of (nb-1)*1,(nb-1)*1,(nb-1)*1,and (nb)*1.
% P,Q,I^2 is arranged in the sequence of P_index; U is arranged in the bus
% index However, note that 点2不能是PV节点.
idgen=zeros(npv,1);
for i = 1: npv
idgen(i,1) =find(gen(:,1)==pv(i,1));
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


% develop U_P矩阵
U_P_var = zeros(nb-1,nb-1+nb-1+nb-1+nb);
for i = 1 :nb-1

       U_P_var(i,nb-1+nb-1+nb-1+P_index(i,1))  = 1;  %根据节点序号来设置电压变量顺序
       U_P_var(i,nb-1+nb-1+nb-1+ P_index(i,2))  = -1;
        
        U_P_var(i,i) = -Lnbr(i,4);
        U_P_var(i,i+nb-1) = -Lnbr(i,5);
        U_P_var(i,i+(nb-1)*2) = Lnbr(i,5)^2 + Lnbr(i,4)^2;

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

delta_P_s = -bus(2:end,3) - [zeros(npq,1);gen(idgen,2)];
delta_Q_s = -bus(2:end,4) ;


% ------------形成Aeq和Beq矩阵-------- %
% Step 3: Form the constant matrix for Aeq*x=beq
Aeq_modify = [ (1/2).*A_P_reduced   zeros(nb-1,nb+npv+nshtc) ;
   (1/2).*A_Q_reduced   zeros(nb-1,nb)   A_UQG A_SHTC ; 
        U_P_var zeros(nb-1,npv+nshtc) ];

Beq_modify = [delta_P_s;delta_Q_s;zeros(nb-1,1)];
% ------------- END ------------------ %
% ------------不等式约束-------------% 
% 电压约束
Aq_U = [zeros(nb,3*(nb-1)) diag(ones(1,nb)) zeros(nb,npv+nshtc)];  
Aq_Qpv = [zeros(npv,3*(nb-1)+nb) diag(ones(npv,1)) zeros(npv,nshtc) ];  
Aq_Qshtc = [zeros(nshtc,3*(nb-1)+nb+npv) diag(ones(nshtc,1))]; 
Aq_QGSW = [zeros(1,nb-1) 1 zeros(1,nb-2 + nb-1+nb+npv+nshtc)];

Aq_max = [Aq_U;Aq_Qpv ;Aq_Qshtc;Aq_QGSW ];

Aq_min = -Aq_max;
Aq = [Aq_max;Aq_min];
bq=[bus(:,9).^2;gen(2:end,5)'; shtc(:,4);gen(1,5)';...
    -bus(:,10).^2;-gen(2:end,4)'; -shtc(:,3);-gen(1,4)']; %线性不等约束f
% ------------- END ------------------ %

% -------------LMI约束（对角化）---------------%
% LMI矩阵约束(SDP约束)
F0 = zeros(2*(nb-1),2*(nb-1)); 

% At_Pij、At_Qij、At_IC和At_U矩阵形成
At_P = cell((nb-1),1);  %nb-1个变量，比如At_P{4}表示Pij第4个变量对应的A4(等于将5个LMI约束组成成的对角化A4矩阵，A4对应第4个变量)
At_Q = cell((nb-1),1);  %nb-1个变量，比如At_P{4}表示Pij第4个变量对应的A4
At_IC = cell((nb-1),1);  %nb-1个变量，比如At_P{4}表示Pij第4个变量对应的A4
At_U = cell(nb,1);  %nb-1个变量，比如At_P{4}表示Pij第4个变量对应的A4
j =sqrt(-1);
for m = 1:nb-1
    
    At_P{m}=zeros(2*(nb-1),2*(nb-1));
    At_Q{m}=zeros(2*(nb-1),2*(nb-1));
    At_IC{m}=zeros(2*(nb-1),2*(nb-1));
    
    for k =1:nb-1
        if m==k
            At_P{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 -1 ;0 0];
            At_Q{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 -j ; 0 0];
            At_IC{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 0 ; 0 -1];
        end
    end
end

for m = 1:nb
    At_U{m}=zeros(2*(nb-1),2*(nb-1));
    for k =1:nb-1
        if m==P_index(k,1)
            At_U{m}((k-1)*2+1:k*2,(k-1)*2+1:k*2) =[-1 0; 0 0];
        end
    end
end


At_P_vec = [];
At_Q_vec = [];
At_IC_vec = [];
for i = 1: nb-1
At_P_vec=[At_P_vec  vec(At_P{i}) ];  
At_Q_vec=[At_Q_vec  vec(At_Q{i}) ];  
At_IC_vec=[At_IC_vec  vec(At_IC{i}) ];  
end
At_U_vec = [];
for i = 1: nb
At_U_vec=[At_U_vec   vec(At_U{i}) ]; 
end
 
% At_Qpv_vec_modify At_shtc_vec_modify
At=[ At_P_vec  At_Q_vec  At_IC_vec At_U_vec  zeros((2.*(nb-1))^2,npv+nshtc)];
Att=[Aeq_modify;-Aq; At ];
ct = [Beq_modify;-bq;vec(F0)];  % Aq*x>bq 写成-bq-（-y*Aq)>0,故前面有负号

% 目标函数
bt=-[1 zeros(1,3*(nb-1)+nb+nshtc+npv-1)]';  % 第1条线路包含平衡节点
 
% ---------SeduMi求解程序----- %
K.f=size(Aeq_modify,1);%等式个数
K.l=size(Aq,1); %非负实数为不等式个数
K.s=size(F0,1);
K.scomplex = 1;
[x_f,y_f,info] = sedumi(Att,bt,ct,K);   
x = y_f;

% verify the LMI constraints
% PP = x(1).*At_P{1}+x(2).*At_P{2}+x(3).*At_P{3} +x(4).*At_P{4}+x(5).*At_P{5} ;
% QQ = x(6).*At_Q{1}+x(7).*At_Q{2}+x(8).*At_Q{3} +x(9).*At_Q{4}+x(10).*At_Q{5} ;
% II = x(11).*At_IC{1}+x(12).*At_IC{2}+x(13).*At_IC{3} +x(14).*At_IC{4}+x(15).*At_IC{5} ;
% UU = x(16).*At_U{1}+x(17).*At_U{2}+x(18).*At_U{3} +x(19).*At_U{4}+x(20).*At_U{5} +x(21).*At_U{6};
% ss = F0 - (PP+QQ+II+UU);
% dd = (ss+ss')/2; %dd is X
% eig(dd)

%-------------------------END ------------------------%
% Nonlinear Constraints（SOC约束）
nlcon = @(x) [ x( (nb-1)*3+P_index(:,1) ).*x( (nb-1)*2+1:(nb-1)*2+nb-1 )  ...
    - ( x( 1: nb-1)./2 ).^2 - ( x(nb-1 + 1:(nb-1)*1 + nb-1)./2 ).^2  ]';     
SOC_value = nlcon(x);

% 输出无功优化结果
busV = sqrt(x(3*(nb-1)+1:3*(nb-1)+nb));  % 最优电压幅值
bus(:,7) = busV;
bus(:,8) = 0;

Pns = x(1,1)./2;                   % 最优注入有功
gen(pv,3) = x(3*(nb-1)+nb+1:3*(nb-1)+nb+npv);
gen(1,3) = x((nb-1)+1)./2;
gen(1,2) = Pns;

shtc(:,2) = x(3*(nb-1)+nb+1,1);           % 最优注入无功补偿
loss = x(1,1)./2+sum(bus(:,3),1);  % 最小网损

Lnbr(:,7) = [x(1:nb-1)./2]';
Lnbr(:,8) = [x((nb-1)+1:(nb-1)+nb-1)./2]';

Lnbr(:,9) = - (x(1:nb-1)./2 - [x((nb-1)*2+1:(nb-1)*2+nb-1)].*Lnbr(:,4));
Lnbr(:,10) = -(x((nb-1)+1:(nb-1)+nb-1)./2 - [x((nb-1)*2+1:(nb-1)*2+nb-1)].*Lnbr(:,5));

Lnbr(:,11) = bus(Lnbr(:,2),7);
Lnbr(:,12) = bus(Lnbr(:,3),7);

I2 =x((nb-1)*2 +1 :(nb-1)*2 + nb-1);

% Step 4: Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if info.feasratio == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

sysdt(PFITER) = info.iter;

% Step 6: Bus Solution and Line Flows 
[gen,bus,Lnbr,trsfm,shtc,shtr,vctr] = int2ext(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss);

elapsed_time = info.cpusec;  %  total optimization time