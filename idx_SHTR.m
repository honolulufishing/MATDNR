function [SHTR_BUS,SHTR_QR,QRMAX,QRMIN,SHTR_SERIES,SHTR_STATE,OPTQR] = idx_SHTR
SHTR_BUS	= 1;  % shunts number
SHTR_QR		= 2;  % reactive power compensate to system (MVAr)(pu)
QRMAX		= 3;  % maximum reactive power (pu)
QRMIN		= 4;  % minimum reactive power (pu)
SHTR_SERIES	= 5;  % shunts series
SHTR_STATE  = 6;  % initial shunt status, 1 - in service, 0 - out of service
OPTQR       = 7;  % optimal reactive power solution for all shunts
return;