function [SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC

SHTC_BUS	= 1;   % 并联电容补偿节点号
SHTC_QC		= 2;   % 并联电容补偿无功(pu)，也是潮流初始值
QCMAX		= 3;   % 无功补偿最大值(pu)
QCMIN		= 4;   % 无功补偿最大值(pu)
SHTC_SERIES	= 5;   % 并联电容组最大组数
SHTC_STATE	= 6;   % 并联电容状态，=1 新建；=0已建
OPTQC		= 7;   % 补偿电容无功功率的优化结果。
return;