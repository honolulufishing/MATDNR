
clc;
clear all;
close all;

% load data
[gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data6_test();
% %[gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data16_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data33;
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data69_test();
% [gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data123_test();
%[gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt] = data1060_test();



% Run ORPF by MOSEK toolbox with SOC relaxation 
Lnbr= Lnbr_all(Lnbr_all(:,7)>=1,:);
[ L_opt,t_opt,SOC_value_opt]= orpf_solver(gen,bus,Lnbr,trsfm,shtc,shtr,vctr,sysdt,'varipm');


L_opt-sum(bus(:,3),1)

t_opt