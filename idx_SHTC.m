function [SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC

SHTC_BUS	= 1;   % �������ݲ����ڵ��
SHTC_QC		= 2;   % �������ݲ����޹�(pu)��Ҳ�ǳ�����ʼֵ
QCMAX		= 3;   % �޹��������ֵ(pu)
QCMIN		= 4;   % �޹��������ֵ(pu)
SHTC_SERIES	= 5;   % �����������������
SHTC_STATE	= 6;   % ��������״̬��=1 �½���=0�ѽ�
OPTQC		= 7;   % ���������޹����ʵ��Ż������
return;