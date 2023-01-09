function [loss,gen,bus_s,Lnbr,trsfm_s,shtc,shtr,vctr,elapsed_time,SOC_value,I2 ] = opti_var_sedumi_new(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt)
% 采用SeduMi与BARON, MOSEK互换模型，即按照MOSEK模型P和Q乘以2
% since 2020.4.02

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

% Step 2:定义优化变量[2*P 2*Q I^2 U^2 QGforPV Qshtc tap^2].其中，P,Q,I^2,U^2 is in the dimensions of (nb-1)*1,(nb-1)*1,(nb-1)*1,and (nb)*1.
% P,Q,I^2 is arranged in the sequence of P_index; U is arranged in the bus
% index However, note that 点2不能是PV节点.
P_index = [frombus tobus];

idgen=zeros(npv,1);
for i = 1: npv
idgen(i,1) =find(gen(:,1)==pv(i,1));
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
U_P_var = zeros(nb-1,nb-1+nb-1+nb-1+nb);
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

delta_P_G = zeros(nb-1,1);
delta_P_G(gen(idgen,1),1) = gen(idgen,2);
delta_P_s = -bus(2:nb,3) - delta_P_G;
delta_Q_s = -bus(2:nb,4) ;

% ------------形成Aeq和Beq矩阵-------- %
% Step 3: Form the constant matrix for Aeq*x=beq
Aeq_modify = [ (1/2).*A_P zeros(nb-1,nb-1) I_P zeros(nb-1,NB_modify+npv+nshtc+ntrsfm) ;
       zeros(nb-1,nb-1) (1/2).*A_Q I_Q zeros(nb-1,NB_modify) A_UQG A_SHTC zeros(nb-1,ntrsfm) ; 
        U_P_var zeros(nb-1,npv+nshtc+ntrsfm)  ];

Beq_modify = [delta_P_s;delta_Q_s;zeros(nb-1,1)];
% ------------- END ------------------ %


% ------------不等式约束-------------% 
% 电压约束
Aq_U = [zeros(NB_modify,3*(nb-1)) diag(ones(1,NB_modify)) zeros(NB_modify,npv+nshtc+ntrsfm)];  
Aq_Qpv = [zeros(npv,3*(nb-1)+NB_modify) diag(ones(npv,1)) zeros(npv,nshtc+ntrsfm) ];  
Aq_Qshtc = [zeros(nshtc,3*(nb-1)+NB_modify+npv) diag(ones(nshtc,1)) zeros(nshtc,ntrsfm)]; 
Aq_Tap = [zeros(nshtc,3*(nb-1)+NB_modify+npv+nshtc)  diag(ones(ntrsfm,1))]; 
Aq_QGSW = [zeros(1,nb-1) 1 zeros(1,nb-2 + nb-1+NB_modify+npv+nshtc+ntrsfm)];

Aq_max = [Aq_U;Aq_Qpv ;Aq_Qshtc;Aq_Tap ;Aq_QGSW ];

Aq_min = -Aq_max;
Aq = [Aq_max;Aq_min];
bq=[0.9.^2*ones(NB_modify,1);gen(2:end,5)'; shtc(:,4);1./trsfm(:,7);2.*gen(1,5)';...
    -1.1.^2*ones(NB_modify,1);-gen(2:end,4)'; -shtc(:,3);-1./trsfm(:,8);-2.*gen(1,4)']; %线性不等约束f

bq(tbnode(:,1),1) = 0.9.^4;
bq(NB_modify+npv+nshtc+ntrsfm+tbnode(:,1),1)=-1.1.^4;
% ------------- END ------------------ %

% -------------LMI约束（对角化）---------------%
% LMI矩阵约束(SDP约束)
F0 = zeros(2*(NB_modify-1),2*(NB_modify-1)); 

% At_Pij、At_Qij、At_IC和At_U矩阵形成
At_P = cell((nb-1),1);  %nb-1个变量，比如At_P{4}表示Pij第4个变量对应的A4(等于将5个LMI约束组成成的对角化A4矩阵，A4对应第4个变量)
At_Q = cell((nb-1),1);  
At_IC = cell((nb-1),1);  
At_U = cell(NB_modify,1);  
At_tap = cell(ntrsfm,1);
j =sqrt(-1);
for m = 1:nb-1
    
    At_P{m}=zeros(2*(NB_modify-1),2*(NB_modify-1));
    At_Q{m}=zeros(2*(NB_modify-1),2*(NB_modify-1));
    At_IC{m}=zeros(2*(NB_modify-1),2*(NB_modify-1));
    
    for k =1:nb-1
        if m==k
            At_P{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 -1 ;0 0];
            At_Q{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 -j ; 0 0];
            At_IC{m}((k-1)*2+1:k*2,(m-1)*2+1:m*2) =[0 0 ; 0 -1];
        end
    end
end

for m = 1  :NB_modify
    At_U{m}=zeros(2*(NB_modify-1),2*(NB_modify-1));
    for k = 1  :nb-1+ntrsfm  %k表示约束序号
        id_tbnode = find(tbnode(:,1)==m); 
        id_tbnode_new = find(tbnode(:,3)==m);
        if k<=nb-1 && (m==P_index(k,1)) && (isempty(id_tbnode))
            At_U{m}((k-1)*2+1:k*2,(k-1)*2+1:k*2) =[-1 0; 0 0];
        elseif  k>nb-1 && (isempty(id_tbnode)) && (~isempty(id_tbnode_new)) 
            At_U{m}((k-1)*2+1:k*2,(k-1)*2+1:k*2) =[0 -2; 0 0];  % tap的LMI约束整合（v7）
        elseif ~isempty(id_tbnode) && k>nb-1 && m==tbnode(id_tbnode,1)
             At_U{m}((k-1)*2+1:k*2,(nb-1+id_tbnode-1)*2+1:(nb-1+id_tbnode)*2) =[-1 0; 0 0];  % tap的LMI约束整合（v1^2）
        end
    end
end

for m = 1:ntrsfm
    At_tap{m}=zeros(2*(NB_modify-1),2*(NB_modify-1));
    for k =nb-1+1:nb-1+ntrsfm
        if nb-1+m==k
            At_tap{m}((k-1)*2+1:k*2,(k-1)*2+1:k*2) =[0 0; 0 -1];
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
At_tap_vec = [];
for i = 1: ntrsfm
    At_tap_vec=[At_tap_vec  vec(At_tap{i}) ];
end
At_U_vec = [];
for i = 1: NB_modify
At_U_vec=[At_U_vec   vec(At_U{i}) ]; 
end

% 组合LMI约束的At
% At_Qpv_vec_modify At_shtc_vec_modify
At=[ At_P_vec  At_Q_vec  At_IC_vec At_U_vec  zeros((2.*(NB_modify-1))^2,npv+nshtc) At_tap_vec ];
Att=[Aeq_modify;-Aq; At ];

ct = [Beq_modify;-bq;vec(F0)];  % Aq*x>bq 写成-bq-（-y*Aq)>0,故前面有负号

% 目标函数
bt=-[1 zeros(1,3*(nb-1)+NB_modify+nshtc+npv+ntrsfm-1)]';  % 第1条线路包含平衡节点
 
% ---------SeduMi求解程序----- %
K.f=size(Aeq_modify,1);%等式个数
K.l=size(Aq,1); %非负实数为不等式个数
K.s=size(F0,1);
K.scomplex = 1;
[x_f,y_f,info] = sedumi(Att,bt,ct,K);   
x = y_f;

% % verify the LMI constraints
% PP = x(1).*At_P{1}+x(2).*At_P{2}+x(3).*At_P{3} +x(4).*At_P{4}+x(5).*At_P{5} ;
% QQ = x(6).*At_Q{1}+x(7).*At_Q{2}+x(8).*At_Q{3} +x(9).*At_Q{4}+x(10).*At_Q{5} ;
% II = x(11).*At_IC{1}+x(12).*At_IC{2}+x(13).*At_IC{3} +x(14).*At_IC{4}+x(15).*At_IC{5} ;
% UU = x(22).*At_U{7}+x(17).*At_U{2}+x(18).*At_U{3} +x(19).*At_U{4}+x(20).*At_U{5} +x(21).*At_U{6};
% ZZ = x(16).*At_U{1}+x(24).*At_tap{1};
% 
% ss = F0 - (PP+QQ+II+UU+ZZ);
% dd = (ss+ss')/2; %dd is X
% eig(dd)

%-------------------------END ------------------------%
% Nonlinear Constraints（SOC约束）
nlcon = @(x) [ x( (nb-1)*3+P_index(:,1) ).*x( (nb-1)*2+1:(nb-1)*2+nb-1 )  ...
    - ( x( 1: nb-1)./2 ).^2 - ( x(nb-1 + 1:(nb-1)*1 + nb-1)./2 ).^2  ]';     
SOC_value = nlcon(x);

% 输出无功优化结果
busV = sqrt(x(3*(nb-1)+1:3*(nb-1)+NB_modify));  % 最优电压幅值
busV(tbnode(:,3),:)=[];  %消去新增节点
bus_s(1:nb,7) = busV;
bus_s(1:nb,8) = 0;
bus_s(tbnode(:,1),7) = sqrt(busV(tbnode(:,1)));

% 输出无功优化结果
busV = sqrt(x(3*(nb-1)+1:3*(nb-1)+NB_modify));  % 最优电压幅值

trsfm_s(:,6) = 1./sqrt(sqrt(x(3*(nb-1)+NB_modify+npv+nshtc+1:3*(nb-1)+NB_modify+npv+nshtc+ntrsfm)));

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

if info.feasratio ~=1
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

sysdt(PFITER) = info.iter;

% Step 6: Bus Solution and Line Flows 
[gen,bus_s,Lnbr,trsfm_s,shtc,shtr,vctr] = int2ext_trs(i2e, gen,bus_s,Lnbr,trsfm_s,shtc,shtr,vctr,sysdt,loss);

elapsed_time = info.cpusec;  %  total optimization time
