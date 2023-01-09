
function [ L_opt,t_opt,SOC_value_opt]= orpf_solver(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,option)
% The DNR Toolbox v1.0 is a Matlab-based package to solve a large-scale DNR problem
% Objective is to Minimize Real Power Loss such that operatioal DistFlow-based constraints
% More can be found at https://github.com/honolulufishing/DNR-MATLAB-Toolbox

% Author(s):Chao Lei
% $Date:2023/01/08 22:19$

if nargin < 9
    error('Please input 9 variables');
end

switch option
    
    case 'socp_mosek'
        % SOCP-based ORPF Model by MOSEK Package
        [loss,gen_s,busV,Lnbr_s,trsfm_s,shtc_s,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_mosek(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt);
        
        L_opt = loss;
        t_opt = elapsed_time;
        SOC_value_opt = SOC_value;
        
       case 'socp_baron'
        % SOCP-based ORPF Model by Baron Package
        [loss,gen_s,busV,Lnbr_s,trsfm_s,shtc_s,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_baron(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt);
        
        L_opt = loss;
        t_opt = elapsed_time;
        SOC_value_opt = SOC_value;   
        
        case 'varipm'
        % Run ORPF for Distribution Networks with VAROPT package
       [bus_s,gen_s,trsfm_s,shtcr_s,loss,elapsed_time,Record_KKT,success] = VARIPM(bus,gen,Lnbr,trsfm,shtc,shtr,vctr);
         
        L_opt = loss;
        t_opt = elapsed_time;
        SOC_value_opt = [];   
        
        case 'sdp_sedumi'
        % SDP-based ORPF Model by SeduMi Package
        [loss,gen_s,busV,Lnbr_s,trsfm_s,shtc_s,shtr,vctr,elapsed_time,SOC_value,I2] = opti_varSDP_sedumi(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt);
       
        L_opt = loss;
        t_opt = elapsed_time;
        SOC_value_opt = SOC_value;   
       
        case 'opti_var_polyhedral_mosek'
        % ORPF Model with Polyhedral Approximation of SOC constraints by MOSEK Package
        [loss,gen_s,busV,Lnbr_s,trsfm_s,shtc_s,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_polyhedral_mosek(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt);
       
        L_opt = loss;
        t_opt = elapsed_time;
        SOC_value_opt = SOC_value;   
end
