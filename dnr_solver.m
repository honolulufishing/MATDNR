
function [ u_opt,L_opt,flag_opt,t_opt,num_branch_opt,SOC_value_opt]= dnr_solver(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt,option)
% The DNR Toolbox v1.0 is a Matlab-based package to solve a large-scale DNR problem
% Objective is to Minimize Real Power Loss such that operatioal DistFlow-based constraints
% More can be found at https://github.com/honolulufishing/DNR-MATLAB-Toolbox

% Author(s):Chao Lei
% $Date:2023/01/08 22:19$

if nargin < 9
    error('Please input 9 variables');
end

switch option
    
    case 'misocp_bigm'
        % MISOCP-based DNR Model with Big-M Relaxation Method by MOSEK Package
        [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc_s,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_BigM(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt = SOC_value;
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end
        
    case 'misocp_mc'
        % MISOCP-based DNR Model using McCormick Linearization Method by MOSEK Package
        [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc_s,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_McCormick(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt = SOC_value;
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end
       
         case 'miqp_bigm'
        % Quadratic DNR Model with Big-M Relaxation Method by MOSEK Package
        [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc_s,shtr,vctr,elapsed_time,num_branch,flag] = opti_DNR_LDF_SCF_BigM_MIQP(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt =[];
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end
        
    case 'milp_bigm'
        % Quadratic DNR Model with Big-M Relaxation Method by MOSEK Package
        [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc_s,shtr,vctr,elapsed_time,num_branch,flag] = opti_DNR_LDF_SCF_BigM(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt =[];
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end
        
    case 'miqp_mc'
        % Quadratic DNR Model with McCormick linearization method Method by MOSEK Package
       [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc,shtr,vctr,elapsed_time,num_branch,flag] = opti_DNR_LDF_SCF_McCormick_MIQP(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
      if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt =[];
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end


        
    case 'polyhedral_bigm'
        % Polyhedral Approximation of DNR Model with Big-M Relaxation Method by MOSEK Package
        [u,loss,gen_s,bus_s,Lnbr_output,trsfm,shtc_s,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_DCH_polyhedron(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt =[];
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end 
        
            case 'misocp_dch'
        % Exact DNR Model with Disjunctive Convex Hull Method by MOSEK Package
        [u,loss,gen,bus,Lnbr_output,trsfm,shtc,shtr,vctr,elapsed_time, SOC_value,num_branch,flag] = opti_DNR_EDF_SCF_DCH(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt);
        if flag==1
            L_opt = loss;
            t_opt = elapsed_time;
            u_opt = u;
            flag_opt = flag;
            num_branch_opt=num_branch;
            SOC_value_opt =[];
        else
            L_opt =  0;
            t_opt = 0;
            u_opt = [];
            flag_opt = 0;
            num_branch_opt=0;
            SOC_value_opt = 0;
        end
end
