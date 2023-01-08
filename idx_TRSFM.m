function [TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,TRS_PF,TRS_QF,TRS_PT,TRS_QT,...
      TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM

TRS_NUM		= 1;        % transformer number
F_BUS		= 2;	    % f, from bus number
T_BUS		= 3;	    % t, to bus number
TRS_R		= 4;	    % r, resistance (p.u.)
TRS_X		= 5;	    % x, reactance (p.u.)
TAP			= 6;	    % ratio, transformer off nominal turns ratio
TAPMAX		= 7;	    % maximum ratio(p.u.)
TAPMIN  	= 8;        % minimum ratio(p.u.)
TAPSERIES	= 9;        % maximum series
TAPSTATE	= 10;       % initial branch status, 1 - in service, 0 - out of service
TRS_PF		= 11;	    % real power injected at "from" bus end (MW)  
TRS_QF		= 12;	    % reactive power injected at "from" bus end (MVAr)
TRS_PT		= 13;	    % real power injected at "to" bus end (MW)  
TRS_QT		= 14;	    % reactive power injected at "to" bus end (MVAr)
TRS_FVM		= 15;       % voltage magnitude at f bus
TRS_TVM		= 16;       % voltage magnitude at t bus
OPTTAP		= 17;       %  
TRS_OPTPF	= 18;	    %
TRS_OPTQF	= 19;	    % 
TRS_OPTPT	= 20;	    % 
TRS_OPTQT	= 21;	    % 
TRS_OPTFVM	= 22;       % 
TRS_OPTTVM	= 23;       %
return;