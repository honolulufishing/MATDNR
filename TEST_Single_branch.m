
clc;
clear all;
close all;

% load data
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data6_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data16_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data33;
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data69_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data123_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data1060_test();
[gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data5_test();

% % Run ORPF by MOSEK toolbox with polyhedral approximation of SOC relaxation 
% Lnbr= Lnbr_all(Lnbr_all(:,7)>=1,:);
% [ L_opt,t_opt,SOC_value_opt]= orpf_solver(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,'socp_mosek');

% Run ORPF by MOSEK toolbox with polyhedral approximation of SOC relaxation 
% Lnbr= Lnbr_all(Lnbr_all(:,7)>=1,:);
% [ L_opt,t_opt,SOC_value_opt]= orpf_solver(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,'socp_mosek');

Lnbr= Lnbr_all(Lnbr_all(:,7)>=1,:);

% Run ORPF for Distribution Networks with VAROPT package
% [bus_s,gen_s,trsfm_s,shtcr_s,loss,elapsed_time,Record_KKT,success] = VARIPM(bus,gen,Lnbr,trsfm,shtc,shtr,vctr);
    
[loss,gen,bus,Lnbr,trsfm,shtc,shtr,vctr,elapsed_time,SOC_value,I2] = opti_var_general_mosek(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt);



