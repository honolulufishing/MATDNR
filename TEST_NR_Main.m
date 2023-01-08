
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

[ u_opt,L_opt,flag_opt,t_opt,num_branch_opt,SOC_value_opt]= dnr_solver(gen,bus,Lnbr_all,trsfm,shtc,shtr,vctr,sysdt,'polyhedral_bigm');

