
function [u,loss,gen,bus,Lnbr_output,trsfm,shtc,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_BigM(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt)
% Exact DistFlow-based Topology Reconfiguration Model using Big-M Method
% Obj: Minimize Real Power Loss over reconfigurable topology and reactive power compensation devices
% s.t. Spanning Tree (ST) Constraints + Single-Commodity Flow (SCF) Constraints + DistFlow equations using Big-M Relxations
% For more information, please refers to our research paper.
% Author: Chao Lei
% e-mail: 21118924r@connnect.polyu.hk
% Date: 2022.08.24

%----------------Pre-Processing-----------------%
% Converts external to internal bus numbering(non-conseutive to consecutive)
[i2e,gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr] = ext2int(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr);

% Intialize and pre-process before using power flow method
[Ybus,Yf,Yt,V0,ref,pv,pq,Sbus,tol,pfmxt] = pfinitial(gen,bus,Lnbr_all,trsfm,shtc,shtr,sysdt);

frombus = Lnbr_all(:,2);
tobus = Lnbr_all(:,3);
nb = size(bus,1); npv=size(pv,1);npq=size(pq,1);
nshtc=size(shtc,1);

nLnbr = size(Lnbr_all,1);

% Optimization vector ids defined as x=[2P 2Q U^2 QG Qshtc f beta u]:
%-------Equality constraints to form matrix Aeq and vector beq--------%
%--------Power balance equality w.r.t P^l,Q^l variables--------%
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
    A_P(i,2*nLnbr+ids) = -2*Lnbr_all(ids,4);
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
    A_Q(i,2*nLnbr+ids) = -2*Lnbr_all(ids,5);
    %  nodes with power flow ejections are negative ones
    idsm =find(P_index(:,1)==i);
    A_Q(i,nLnbr+idsm) = -1;
end
A_Q_reduced = A_Q(2:nb,1:3*nLnbr);

% construct the coincidence QG matrix
A_UQG = zeros(nb-1,npv);
for i =1:npv
    A_UQG(pv(i,1)-1,i) =1;
end

% construct the coincidence shtc matrix
A_SHTC = zeros(nb-1,nshtc);
for i = 1 :nshtc
    A_SHTC(shtc(i,1)-1,i) =1;
end

delta_P_s = -bus(2:end,3) - [zeros(npq,1);gen(idgen,2)];
delta_Q_s = -bus(2:end,4) ;

% Add auxiliary variables for SOC Constraints
% auxiliary variable w for equality constraints
A_auxiliary_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1 :nLnbr
    A_auxiliary_1(i,2*nLnbr+i) = 1;
    A_auxiliary_1(i,2*nLnbr+nLnbr+P_index(i,1)) =1;
    A_auxiliary_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+i) = -1;
end

% auxiliary variable m for equality constraints
A_auxiliary_2 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1 : nLnbr
    A_auxiliary_2(i,2*nLnbr+i) = 1;
    A_auxiliary_2(i,2*nLnbr+nLnbr+P_index(i,1)) =-1;
    A_auxiliary_2(i,2*nLnbr+nLnbr+nb+npv+nshtc+nLnbr+i) = -1;
end
% ------------- END ------------------ %

% --------Radiality Constraints w.r.t beta variables---------%
A_f = zeros(nb,nLnbr);
for i = 2 : nb
    % nodes with power flow injected to A_P are ones
    ids =find(P_index(:,2)==i);
    A_f(i,ids) = 1;
    %  nodes with power flow ejected out of A_P are negative ones
    idsm =find(P_index(:,1)==i);
    A_f(i,idsm) = -1;
end
b_f = [0;ones(nb-1,1)];

% inequality constraints for u^l*Smin<f^l< u^l*Smax
S_max = nb;
A_f_u_1 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
A_f_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1: nLnbr
    A_f_u_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+i) = 1;
    A_f_u_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+i) = 1;
end
A_f_u_1(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=-S_max*diag(ones(1,nLnbr));
A_f_u_2(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=S_max*diag(ones(1,nLnbr));

A_beta_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1:nLnbr
    A_beta_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+(i-1)*2+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+2*i) =[1 1];
    A_beta_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+i) =-1;
end
b_beta_1 = zeros(nLnbr,1);

% nodal constraints for source nodes
n_source = size(find(P_index(:,1)==1),1) + size(find(P_index(:,2)==1),1);
A_beta_2 = zeros(n_source,2*nLnbr);
idf = find(P_index(:,1)==1);   % 1 is the reference node
for i = 1 :  n_source
    A_beta_2(i,(idf(i,1)-1).*2+1) =1;
end
b_beta_2 = zeros(n_source ,1);

% For directional flows, switch status can be explicitly accommodated
n_ex = size(Lnbr_all(Lnbr_all(:,7)==2,1),1);
A_beta_3 = zeros(2*n_ex,2*nLnbr);
b_beta_3 = zeros(2*n_ex,1);
ide = find(Lnbr_all(:,7)==2);
for i = 1 :  n_ex
    A_beta_3((i-1)*2+1,(ide(i,1)-1).*2+1) = 1;
    b_beta_3((i-1)*2+1,1) = 0;
    
    A_beta_3(2*i,ide(i,1)*2) = 1;
    b_beta_3(2*i,1) = 1;
end

% nodal constraints for PQ loads
A_beta_4 = zeros(nb-1,2*nLnbr);
b_beta_4 = ones(nb-1,1);
for i = 2 : nb
    idl = find(P_index(:,1)==i);
    for j = 1 : size( idl ,1)
        if Lnbr_all(idl(j,1),7)~=2
        A_beta_4(i-1,idl(j,1)*2-1)=1;  
         end
    end
    
    idt = find(P_index(:,2)==i);
    for j = 1 : size( idt ,1)
        A_beta_4(i-1,idt(j,1)*2)=1;
    end
end

A_radiality = ones(1,nLnbr);
b_radiality = nb-1;

A_beta = [zeros(nb,3*nLnbr+nb+npv+nshtc+2*nLnbr) A_f zeros(nb,3*nLnbr);
    A_beta_1 ;
    zeros(n_source,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr) A_beta_2 zeros(n_source,nLnbr);
    zeros(2*n_ex,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr) A_beta_3 zeros(2*n_ex,nLnbr);
    zeros(nb-1,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr) A_beta_4 zeros(nb-1,nLnbr);
   zeros(1,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr) A_radiality  
  ];

B_beta = [
    b_f;
    b_beta_1 ; 
    b_beta_2 ;
    b_beta_3; 
    b_beta_4 ;
       b_radiality;
    ];
% -------------END------------------ %

%----Form the constant matrix for Aeq*x=beq-----%
Aeq_modify = [ (1/2).*A_P_reduced       zeros(nb-1,nb+npv+nshtc+2*nLnbr)       zeros(nb-1,nLnbr*4);
    (1/2).*A_Q_reduced       zeros(nb-1,nb)    A_UQG A_SHTC  zeros(nb-1,2*nLnbr+nLnbr*4);
    A_beta ;
    A_auxiliary_1;
    A_auxiliary_2];

Beq_modify = [delta_P_s;delta_Q_s;B_beta;zeros(2*nLnbr,1)];
% ------------- END ------------------ %

%----Specify the cones-----%
prob.cones   = cell(nLnbr,1);
for i = 1 : nLnbr
    prob.cones{i}.type = 'MSK_CT_QUAD';%MSK_CT_RQUAD for Qudratic Cone
    prob.cones{i}.sub  = [3*nLnbr+nb+nshtc+npv+i,  i,nLnbr+i,3*nLnbr+nb+nshtc+npv+nLnbr+i];
end
% ------------- END ------------------ %

% ---------Inequality constraints to form matrix A and vector b------------%
% inequality constraints for -Uf+Ut+2PR+2QX-I*z^2=0
M = 1000;
U_P_var_1 = zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
b_P_var_1 = zeros(nLnbr,1);
for i = 1 :nLnbr
    if Lnbr_all(i,7) == 2
        U_P_var_1(i,3*nLnbr+ frombus(i,1))  = -1;  
        U_P_var_1(i,3*nLnbr+ tobus(i,1))  = 1;
        
        U_P_var_1(i,i) = Lnbr_all(i,4);  % equal to 2PR and note that cause Pl is positive
        U_P_var_1(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_P_var_1(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
    end
    
    if Lnbr_all(i,7)~= 2
        U_P_var_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+i) = M;
        b_P_var_1(i,1)=M;
        
        U_P_var_1(i,3*nLnbr+ frombus(i,1))  = -1; 
        U_P_var_1(i,3*nLnbr+ tobus(i,1))  = 1;
        
        U_P_var_1(i,i) = Lnbr_all(i,4);  % equal to 2PR and note that cause Pl is positive
        U_P_var_1(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_P_var_1(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
    end
end

U_P_var_2 = zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
b_P_var_2= zeros(nLnbr,1);
for i = 1 :nLnbr
    if Lnbr_all(i,7) == 2
        U_P_var_2(i,3*nLnbr+ frombus(i,1))  = -1;  %根据节点序号来设置电压变量顺序
        U_P_var_2(i,3*nLnbr+ tobus(i,1))  = 1;
        
        U_P_var_2(i,i) = Lnbr_all(i,4);
        U_P_var_2(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_P_var_2(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
    end
    
    if Lnbr_all(i,7)~= 2
        U_P_var_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+i) = -M;
        b_P_var_2(i,1)=-M;
        
        U_P_var_2(i,3*nLnbr+ frombus(i,1))  = -1;  %根据节点序号来设置电压变量顺序
        U_P_var_2(i,3*nLnbr+ tobus(i,1))  = 1;
        
        U_P_var_2(i,i) = Lnbr_all(i,4);
        U_P_var_2(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_P_var_2(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
    end
end

% inequality constraints for u^l*Smin<P^l,Q^l,I^l< u^l*Smax
S_max = 1000;
A_P_u_1 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
A_P_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1: nLnbr
    A_P_u_1(i,i) = 1;
    A_P_u_2(i,i) = 1;
end
A_P_u_1(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+2*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=-S_max*diag(ones(1,nLnbr));
A_P_u_2(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+2*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=S_max*diag(ones(1,nLnbr));

A_Q_u_1 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
A_Q_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1: nLnbr
    A_Q_u_1(i,nLnbr+i) = 1;
    A_Q_u_2(i,nLnbr+i) = 1;
end
A_Q_u_1(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+2*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=-S_max*diag(ones(1,nLnbr));
A_Q_u_2(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+2*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=S_max*diag(ones(1,nLnbr));


A_I_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr);
for i = 1: nLnbr
    A_I_u_2(i,2*nLnbr+i) = 1;
end
A_I_u_2(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr)=S_max*diag(ones(1,nLnbr));

% Linear Constraints of QG 
% Reactive power capacity constraints for the reference bus 1
A_QGSW_inequal =  [zeros(1,nLnbr) 1 zeros(1,nLnbr-1) zeros(1,nLnbr) zeros(1,nb) zeros(1,npv) zeros(1,nshtc) zeros(1,2*nLnbr) zeros(1,4*nLnbr)];
b_QGSW_inequal_upper = gen(1,4);
b_QGSW_inequal_lower = gen(1,5);
% ------------- END ------------------ %

% --------------Solver: MOSEK-------------%
% Solve with MOSEK toolbox
% Bounds
lb = [-1000.*ones(1,3*nLnbr) (bus(:,9).^2)' gen(2:end,5)'  shtc(:,4)' -1000.*ones(1,2*nLnbr)   -nb*2.*ones(1,nLnbr) zeros(1,2*nLnbr), zeros(1,nLnbr) ];
ub =[1000.*ones(1,3*nLnbr) (bus(:,10).^2)' gen(2:end,4)'  shtc(:,3)'  1000.*ones(1,2*nLnbr)      nb*2.*ones(1,nLnbr) ones(1,2*nLnbr) ones(1,nLnbr) ];

A = [ Aeq_modify;
    U_P_var_1;
    U_P_var_2;
    A_QGSW_inequal;
    A_P_u_1;
    A_P_u_2;
    A_Q_u_1;
    A_Q_u_2;
    A_I_u_2;
    A_f_u_1;
    A_f_u_2;
    ];
rl = [ Beq_modify;
    -2000.*ones(nLnbr,1);
    b_P_var_2;
    b_QGSW_inequal_lower;
    -2000.*ones(nLnbr,1);
    zeros(nLnbr,1);
    -2000.*ones(nLnbr,1);
    zeros(nLnbr,1);
    zeros(nLnbr,1);
    -2*nb.*ones(nLnbr,1);
    zeros(nLnbr,1);];
ru = [Beq_modify ;
    b_P_var_1;
    2000.*ones(nLnbr,1);
    b_QGSW_inequal_upper;
    zeros(nLnbr,1);
    2000.*ones(nLnbr,1);
    zeros(nLnbr,1);
    2000.*ones(nLnbr,1);
    2000.*ones(nLnbr,1);
    zeros(nLnbr,1);
    2*nb.*ones(nLnbr,1);];

% Specify the non-conic part of the problem.
id_source = find(P_index(:,1)==1);
obj_1 = zeros(1,nLnbr);
obj_1(1,id_source)=1;

prob.c = [obj_1 zeros(1,nLnbr*2+nb+nshtc+npv+2*nLnbr+4*nLnbr) ];

prob.a = sparse(A);
prob.blc = sparse(rl);
prob.buc = sparse(ru);
prob.blx = lb;
prob.bux = ub;
prob.ints.sub=[3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+4*nLnbr];   % u is a set of 0-1 integer variable

% Optimize the problem.
[r,res]=mosekopt('minimize info',prob);
%-------------------------END ------------------------%

% --------------Output: Optimal Solution-------------%
% Display the primal solution.
X_Opt = res.sol.int.xx;

% Output: Optimal bus voltage profile
busV = sqrt(X_Opt(3*nLnbr+1:3*nLnbr+nb));
bus(:,7) = busV;
bus(:,8) = 0;

% Output: Minimal real power injection at the root bus
Pns = X_Opt(1,1)./2;
gen(pv,3) = X_Opt(3*nLnbr+nb+1:3*nLnbr+nb+npv);
gen(1,3) = X_Opt(nLnbr+1)./2;
gen(1,2) = Pns;
gen(1,6) = sqrt(X_Opt(3*nLnbr+1,1));

% Output: Optimal reactive power compensation
shtc(1:nshtc,2) = X_Opt(3*nLnbr+nb+npv+1:3*nLnbr+nb+npv+nshtc,1);
% Output: Minimal real power loss
loss = X_Opt(1,1)./2+sum(bus(:,3),1);

Lnbr(:,1:6) = Lnbr_all(:,1:6);
Lnbr(:,7) = [X_Opt(1:nLnbr)./2]';
Lnbr(:,8) = [X_Opt(nLnbr+1:nLnbr+nLnbr)./2]';

Lnbr(:,9) = - X_Opt(1:nLnbr)./2 ;
Lnbr(:,10) = -X_Opt(nLnbr+1:nLnbr+nLnbr)./2;

Lnbr(:,11) = bus(Lnbr_all(:,2),7);
Lnbr(:,12) = bus(Lnbr_all(:,3),7);

id_status = find( Lnbr_all(:,7)==2);

% Output: Optimal switch status u^l of each branch
u = X_Opt(3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr+1:3*nLnbr+nb+nshtc+npv+2*nLnbr+4*nLnbr,1);   % u is a set of 0-1 integer variable
beta = X_Opt(3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+3*nLnbr,1);   % u is a set of 0-1 integer variable

Lnbr(:,13) = u;
Lnbr(id_status,13) = 2;

% Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if res.rcode == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

% Nonlinear Constraints（SOC Constraints）
nlcon = @(x) [ x( nLnbr*3+Lnbr(:,2)).*x( nLnbr*2+1:nLnbr*3)  ...
    - ( x( 1: nLnbr)./2 ).^2 - ( x(nLnbr + 1:nLnbr + nLnbr)./2 ).^2  ]';
SOC_value = max(nlcon(X_Opt)');
% SOC =  X_Opt( nLnbr*3+Lnbr(:,2)).*X_Opt(nLnbr*2+1:nLnbr*3)  - ( X_Opt( 1: nLnbr)./2 ).^2 - ( X_Opt(nLnbr + 1:nLnbr + nLnbr)./2 ).^2

% output solution
% PF with newton-raphson method
[gen,bus_exact,Lnbr_s,trsfm,shtc,shtr,vctr,loss_exact] = pf(gen,bus,Lnbr(Lnbr(:,13)>0,1:6),trsfm,shtc,shtr,vctr,sysdt);
if max(abs( bus_exact(:,7) - bus(:,7)  ))<0.01
    flag = 1;
    loss = loss_exact;
else
    flag = 0;
end

% Bus Solution and Line Flows
[gen,bus,Lnbr,trsfm,shtc,shtr,vctr] =  int2ext_NR(i2e, gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,loss,beta);

Lnbr_output = Lnbr;
elapsed_time = res.info.MSK_DINF_OPTIMIZER_TIME;  %  total optimization time
num_branch = res.info.MSK_IINF_MIO_NUM_BRANCH;  %  total branches of B&B algorithm
%-------------------------END ------------------------%