function [u,loss,gen,bus,Lnbr_output,trsfm,shtc,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_DCH(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt)
% Exact DistFlow-based Topology Reconfiguration Model using Disjunctive Convex Hull
% Note: Conventinal Convex Hull for (1c)
% Obj: Minimize Real Power Loss over reconfigurable topology and reactive power compensation devices
% s.t. Spanning Tree (ST) Constraints + Single-Commodity Flow (SCF) Constraints + DistFlow equations using disjunctive convex hull relaxations
% where quadratic equations in DistFlow equations are tightned by convex
% hull relaxations for a single branch.
% For more information, please refers to our research paper.
% Author: Chao Lei
% e-mail: 21118924r@connnect.polyu.hk
% Date: 2022.08.21


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
A_auxiliary_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1 :nLnbr
    A_auxiliary_1(i,2*nLnbr+i) = 1;
    A_auxiliary_1(i,2*nLnbr+nLnbr + P_index(i,1)) =1;
    A_auxiliary_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+i) = -1;
end

% auxiliary variable m for equality constraints
A_auxiliary_2 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1 : nLnbr
    A_auxiliary_2(i,2*nLnbr+i) = 1;
    A_auxiliary_2(i,2*nLnbr+nLnbr+ P_index(i,1)) =-1;
    A_auxiliary_2(i,2*nLnbr+nLnbr+nb+npv+nshtc+nLnbr+i) = -1;
end

A_I = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1 : nLnbr
    A_I(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+i) = 1;
    A_I(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+i) =1;
    A_I(i,2*nLnbr+i) = -1;
end

% equality constraints for y1+(2PR+2QX-I*z^2)=0
U_var = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
b_var = zeros(nLnbr,1);
for i = 1 :nLnbr
    
    if Lnbr_all(i,7) == 2
        U_var(i,3*nLnbr+ frombus(i,1))  = -1;
        U_var(i,3*nLnbr+ tobus(i,1))  = 1;
        
        U_var(i,i) = Lnbr_all(i,4);  % equal to 2PR and note that cause Pl is positive
        U_var(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_var(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
    end
    
    if Lnbr_all(i,7)~= 2
        U_var(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i)  =-1;  % y1
        
        U_var(i,i) = Lnbr_all(i,4);  % equal to 2PR and note that cause Pl is positive
        U_var(i,i+nLnbr) = Lnbr_all(i,5);
        
        U_var(i,2*nLnbr+i) = -(Lnbr_all(i,4).^2+Lnbr_all(i,5).^2);  % equal to I*z^2
        
    end
end

% alpha equality for (alpha+)+(alpha_)=beta_nm
alpha_var_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
b_alpha_var_1 = zeros(nLnbr,1);
for i =1 :nLnbr
    
    if Lnbr_all(i,7)~= 2
        alpha_var_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+i)  =1;  % alpha+
        alpha_var_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+nLnbr+i)  =1;  % alpha-
        alpha_var_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+i)  =-1;  % u^l
        b_alpha_var_1(i,1)=0;
    end
    
end


% id=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,33,34,35,36];
id_status = find( Lnbr_all(:,7)==2);
id_remove =[];
for i = 1 : nshtc
    
    id_f = find(P_index(:,1)==shtc(i,1));
    id_t = find(P_index(:,2)==shtc(i,1));
    
    if ~isempty(id_f)
        id_remove =  [id_remove ; id_f];
    elseif ~isempty(id_t)
        id_remove =  [id_remove ; id_t];
    end 
end

id_remove_2 =[];
for i = 1 : size(id_remove,1)
    
    id_f = find(P_index(:,2)==P_index(id_remove(i,1),1));
    id_t = find(P_index(:,2)==P_index(id_remove(i,1),2));
    
    if ~isempty(id_f)
        id_remove_2 =  [id_remove_2 ; id_f];
    elseif ~isempty(id_t)
        id_remove_2 =  [id_remove_2 ; id_t];
    end 
end
id = setdiff(Lnbr_all(:,1),[id_status;id_remove;id_remove_2]);


alpha_beta_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+4*nLnbr);
for k = 1: length(id)
    i = id(k);
    alpha_beta_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+i) = 1;
    alpha_beta_1(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*i) = -1;
end

alpha_beta_2 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+4*nLnbr);
for k = 1: length(id)
     i = id(k);
    alpha_beta_2(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+nLnbr+i) = 1;
    alpha_beta_2(i,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*i-1) = -1;
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
A_f_u_1 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
A_f_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1: nLnbr
    A_f_u_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = 1;
    A_f_u_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = 1;
end
A_f_u_1(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr)=-S_max*diag(ones(1,nLnbr));
A_f_u_2(:,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr)=S_max*diag(ones(1,nLnbr));

A_beta_1 = zeros(nLnbr,2*nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1:nLnbr
    A_beta_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+(i-1)*2+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*i) =[1 1];
    A_beta_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+i) =-1;
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

% radial constraints for branches at a number of nb-1
A_radiality = ones(1,nLnbr);
b_radiality = nb-1;

A_beta = [zeros(nb,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr) A_f zeros(nb,3*nLnbr+2*nLnbr+2*nLnbr);
    A_beta_1 ;
    zeros(n_source,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr) A_beta_2 zeros(n_source,nLnbr+2*nLnbr+2*nLnbr);
    zeros(2*n_ex,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr) A_beta_3 zeros(2*n_ex,nLnbr+2*nLnbr+2*nLnbr);
    zeros(nb-1,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr) A_beta_4 zeros(nb-1,nLnbr+2*nLnbr+2*nLnbr);
    zeros(1,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*nLnbr) A_radiality zeros(1,2*nLnbr+2*nLnbr);
    ];

B_beta = [b_f;
    b_beta_1 ;
    b_beta_2 ;
    b_beta_3;
    b_beta_4 ;
    b_radiality;
    ];
% -------------END------------------ %

%----Form the constant matrix for Aeq*x=beq-----%
Aeq_modify = [ (1/2).*A_P_reduced       zeros(nb-1,nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr)       zeros(nb-1,nLnbr*4+2*nLnbr+2*nLnbr);
    (1/2).*A_Q_reduced       zeros(nb-1,nb)    A_UQG A_SHTC  zeros(nb-1,2*nLnbr+2*nLnbr+nLnbr+nLnbr*4+2*nLnbr+2*nLnbr);
    A_beta ;
    A_auxiliary_1;
    A_auxiliary_2;
    A_I;
    U_var ;
    alpha_var_1;
       alpha_beta_1;
       alpha_beta_2;
    ];

Beq_modify = [delta_P_s;
    delta_Q_s;
    B_beta;
    zeros(2*nLnbr,1);
    zeros(nLnbr,1);
    b_var;
    b_alpha_var_1;
       zeros(2*nLnbr,1);
    ];
% ------------- END ------------------ %


%----Specify the cones-----%
prob.cones   = cell(nLnbr,1);
for i = 1 : nLnbr
    prob.cones{i}.type = 'MSK_CT_QUAD';%MSK_CT_RQUAD for Qudratic Cone
    prob.cones{i}.sub  = [3*nLnbr+nb+nshtc+npv+i,  i,nLnbr+i,3*nLnbr+nb+nshtc+npv+nLnbr+i];
end
% ------------- END ------------------ %

% ---------Inequality constraints to form matrix A and vector b------------%
%-------Relaxation constraints for y1-----------%
A_y_1 = zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
A_y_2 = zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
b_y_1 = zeros(nLnbr,1);b_y_2 = zeros(nLnbr,1);
for i = 1 : nLnbr
    if Lnbr_all(i,7)~= 2
        % y1 >= ( v_m - v_n) + (u^l-1)*(1.1^2-0.9^2)
        A_y_1(i,3*nLnbr+ P_index(i,1))  = -1;
        A_y_1(i,3*nLnbr+P_index(i,2))  = 1;
        A_y_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = 1;
        A_y_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+i) =-( bus(Lnbr_all(i,2),10).^2 - bus(Lnbr_all(i,3),9).^2);
        b_y_1(i,1) = -(bus(Lnbr_all(i,2),10).^2 - bus(Lnbr_all(i,3),9).^2);
        
        % y1 <= (v_m - v_n) + (u^l-1)*(0.9^2-1.1^2)
        A_y_2(i,3*nLnbr+P_index(i,1))  = -1;
        A_y_2(i,3*nLnbr+P_index(i,2))  = 1;
        A_y_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = 1;
        A_y_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+i) = -(bus(Lnbr_all(i,2),9).^2 - bus(Lnbr_all(i,3),10).^2);
        b_y_2(i,1) = -(bus(Lnbr_all(i,2),9).^2 - bus(Lnbr_all(i,3),10).^2);
        
    end
end
% ------End-----------%

%---------Perspective cuts for I-----%
S_max=100;
A_I_u_1 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
A_I_u_2 =zeros(nLnbr,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
for i = 1: nLnbr
    A_I_u_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+i) = 1;
    A_I_u_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*i)=-S_max;
    
    A_I_u_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+i) = 1;
    A_I_u_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+2*i-1)=-S_max;
end
% ------End-----------%

%---------Disjunctive cuts for v and y1-----%
A_P_1 = zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
b_P_1 = zeros(nLnbr,1);
A_P_2 = zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
b_P_2 = zeros(nLnbr,1);

A_P_5 = zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
A_P_6 = zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);

A_P_7 = zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);
A_P_8= zeros(nLnbr,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr);

for k = 1: length(id)
    i = id(k);
    if Lnbr_all(i,7)~= 2
        
        % beta_nm: v_m - v_n >= (1 - alpha_+)* delta_v_min
        A_P_1(i,3*nLnbr + P_index(i,1)) = 1;
        A_P_1(i,3*nLnbr + P_index(i,2)) = -1;
        A_P_1(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+i) = -(1.1^2-0.9^2);
        b_P_1(i,1) = -(1.1^2-0.9^2);
        
        % beta_nm: v_m - v_n <= (1 - alpha_-)* delta_v_max
        A_P_2(i,3*nLnbr + P_index(i,1)) = 1;
        A_P_2(i,3*nLnbr + P_index(i,2)) = -1;
        A_P_2(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+nLnbr+i) = (1.1^2-0.9^2);
        b_P_2(i,1) = (1.1^2-0.9^2);
       
        % perspective cuts for y1
        % y1 <= alpha_+*delta_v_max
        A_P_5(i,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = 1;
        A_P_5(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+i) = -(1.1^2-0.9^2);
        % -y1 <= alpha_-*delta_v_max
        A_P_6(i,nLnbr+nLnbr+nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+i) = -1;
        A_P_6(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+nLnbr+i) = -(1.1^2-0.9^2);

        % valid cuts for SOC constraints
        A_P_7(i,i) = Lnbr_all(i,4);
        A_P_7(i,nLnbr+i) = Lnbr_all(i,5);
        A_P_7(i,2*nLnbr+i) = -(Lnbr_all(i,5).^2 + Lnbr_all(i,4).^2);
        A_P_7(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+i) = -(1.1^2-0.9^2);
        
        A_P_8(i,i) = -Lnbr_all(i,4);
        A_P_8(i,nLnbr+i) = -Lnbr_all(i,5);
        A_P_8(i,2*nLnbr+i) = (Lnbr_all(i,5).^2 + Lnbr_all(i,4).^2);
        A_P_8(i,3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+nLnbr+i) = -(1.1^2-0.9^2);
        
    end
end
% -------------------End------------------%

% Linear Constraints of QG ranged between lower and upper boundaries
A_QGSW_inequal =  [zeros(1,nLnbr) 1 zeros(1,nLnbr-1) zeros(1,nLnbr) zeros(1,nb) zeros(1,npv) zeros(1,nshtc) zeros(1,2*nLnbr+2*nLnbr+nLnbr) zeros(1,4*nLnbr+2*nLnbr+2*nLnbr)];
b_QGSW_inequal_upper = gen(1,4);
b_QGSW_inequal_lower = gen(1,5);
% ------------- END ------------------ %

% --------------Solver: MOSEK-------------%
% Solve with MOSEK toolbox
% Bounds
lb = [-100.*ones(1,2*nLnbr) zeros(1,nLnbr) (bus(:,9).^2)' gen(2:end,5)'  shtc(:,4)' -100.*ones(1,2*nLnbr)  zeros(1,2*nLnbr)   -(1.1^2-0.9^2).*ones(1,nLnbr)   -nb*2.*ones(1,nLnbr)  zeros(1,2*nLnbr), zeros(1,nLnbr+2*nLnbr+2*nLnbr) ];
ub =[100.*ones(1,3*nLnbr) (bus(:,10).^2)' gen(2:end,4)'  shtc(:,3)'  100.*ones(1,2*nLnbr) nb.*ones(1,2*nLnbr)  (1.1^2-0.9^2).*ones(1,nLnbr) nb*2.*ones(1,nLnbr) ones(1,2*nLnbr) ones(1,nLnbr+2*nLnbr+2*nLnbr) ];

A = [ Aeq_modify;
    A_QGSW_inequal;
    A_f_u_1;
    A_f_u_2;
    A_I_u_1;
    A_I_u_2;
    A_y_1;
    A_y_2;
    
    A_P_1;
    A_P_2;

    A_P_5;
    A_P_6;
    A_P_7;
    A_P_8;
    ];

rl = [ Beq_modify;
    b_QGSW_inequal_lower;
    -2*nb.*ones(nLnbr,1);
    zeros(nLnbr,1);
    -200.*ones(nLnbr,1);
    -200.*ones(nLnbr,1);
    b_y_1 ;
    -200.*ones(nLnbr,1);
    
    b_P_1 ;
    -200.*ones(nLnbr,1);

    -200.*ones(nLnbr,1);
    -200.*ones(nLnbr,1);
    -200.*ones(nLnbr,1);
    -200.*ones(nLnbr,1);
    
    ];

ru = [Beq_modify ;
    b_QGSW_inequal_upper;
    zeros(nLnbr,1);
    2*nb.*ones(nLnbr,1);
    zeros(nLnbr,1);
    zeros(nLnbr,1);
    200.*ones(nLnbr,1);
    b_y_2;
    
    200.*ones(nLnbr,1);
    b_P_2 ;

    zeros(nLnbr,1);
    zeros(nLnbr,1);
    zeros(nLnbr,1);
    zeros(nLnbr,1);
    ];

% If this option is turned on outer approximation is used when solving relaxations of conic problems; otherwise interior point is used.
% param.MSK_IPAR_MIO_CONIC_OUTER_APPROXIMATION = 1;

% Specify the non-conic part of the problem.
id_source = find(P_index(:,1)==1);
obj_1 = zeros(1,nLnbr);
obj_1(1,id_source)=1;

prob.c = [obj_1 zeros(1,nLnbr*2+nb+nshtc+npv+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr+2*nLnbr+2*nLnbr) ];
prob.a = sparse(A);
prob.blc = sparse(rl);
prob.buc = sparse(ru);
prob.blx = lb;
prob.bux = ub;
prob.ints.sub=[3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+2*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr];   % u is a set of integer variables

% Specify start guess for the integer variables.
% prob.sol.int.xx = [nan.*ones(3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+3*nLnbr,1);Lnbr_all(:,7)]';

% Request constructing the solution from integer variable values
% If set to MSK ON and all integer variables have been given a value for which a feasible mixed
% integer solution exists, then MOSEK generates an initial solution to the mixed integer
% problem by fixing all integer values and solving the remaining problem
% param.MSK_IPAR_MIO_CONSTRUCT_SOL = 1;

[r,res] = mosekopt('minimize info',prob);
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


% Output: Optimal switch status u^l of each branch
u = X_Opt(3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+4*nLnbr,1);   % u is a set of 0-1 integer variable
beta = X_Opt(3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+2*nLnbr+nLnbr+3*nLnbr,1);   % u is a set of 0-1 integer variable


Lnbr(:,13) = u;
Lnbr(id_status,13) = 2;

% Record the Solver Performance Information
[DPRATE,INTRATE,OMRATE,AUECOST,AUCCOST,AURCOST,TMAX,BASEMVA,PFMETHOD,OPTMODEL,OPTMETHOD,...
    ACCURACY,PFMAXIT,OPFMAXIT,POPNUM,CPOPT,TARGET,SUCCESS,PFITER,OPFITER] = idx_sysdt;

if res.rcode == 0
    sysdt(SUCCESS) = 0;
else sysdt(SUCCESS) = 1;
end

% Nonlinear Constraints??¡ìSOC Constraints??
nlcon = @(x) [ x( nLnbr*3+Lnbr(:,2)).*x( nLnbr*2+1:nLnbr*3)  ...
    - ( x( 1: nLnbr)./2 ).^2 - ( x(nLnbr + 1:nLnbr + nLnbr)./2 ).^2  ]';
SOC_value = max(nlcon(X_Opt)');

% output solution
% PF with newton-raphson method
[gen,bus_exact,Lnbr_s,trsfm,shtc,shtr,vctr,loss_exact] = pf(gen,bus,Lnbr(u>0,1:6),trsfm,shtc,shtr,vctr,sysdt);
if max(abs( bus_exact(:,7)-bus(:,7)  ))<0.01
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

% X_l = X_Opt(1:2*nLnbr,1)./2;
% X_u = X_Opt(3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+nLnbr+1:3*nLnbr+nb+npv+nshtc+2*nLnbr+nLnbr+nLnbr+2*nLnbr,1);
%
% plot_convex_hull_EDF(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt,X_l,X_u)
