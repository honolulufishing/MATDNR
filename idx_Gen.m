function [GEN_BUS,GEN_PG,GEN_QG,QGMAX,QGMIN,GEN_VG,GEN_OPTVG,GEN_OPTQG,GEN_OPTPG] = idx_Gen
%IDX_BRCH   Defines variables for column indices to Gen.
GEN_BUS	   = 1;	%% bus number
GEN_PG	   = 2;	%% Pg, real power output (pu)
GEN_QG	   = 3;	%% Qg, reactive power output (pu)
QGMAX	   = 4;	%% Qmax, maximum reactive power output (pu)
QGMIN	   = 5;	%% Qmin, minimum reactive power output (pu)
GEN_VG	   = 6;	%% Vg, voltage magnitude setpoint (pu)
GEN_OPTVG  = 7;   %  the result of after reactive optimization
GEN_OPTQG  = 8;   %  the result of after reactive optimization
GEN_OPTPG  = 9;  %  the result of after reactive optimization
return;
