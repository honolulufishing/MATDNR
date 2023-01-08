function [BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus
%IDX_BRCH   Defines variables for column indices to Bus.
BUS_I	  	= 1;	% bus number (1 to 29997)
BUS_TYPE	= 2;	% bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus
PD			= 3;	% real power demand (pu)
QD			= 4;	% reactive power demand (pu)
GS			= 5;	% shunt conductance (pu)
BS			= 6;	% shunt susceptance (p.u.)
VM			= 7;	% voltage magnitude (p.u.),the result of power flow
VA			= 8;	% voltage angle (degrees),the result of power flow
BASE_KV	    = 9;	% base voltage (kV)
VMAX		= 10;	% maximum voltage magnitude (p.u.)	
VMIN		= 11;	% minimum voltage magnitude (p.u.)	
%% included in power flow solution, not necessarily in input
OPTVM    = 12; %  voltage magnitude after reactive optimization
OPTVA    = 13; %  voltage angle after reactive optimization
return;
