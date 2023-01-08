function [LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,LN_OPTPF,LN_OPTQF,...
      LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR
LN_NUM		    = 1;    % line branches number
F_BUS			= 2;	% f, from bus number
T_BUS			= 3;	% t, to bus number
LN_R			= 4;	% r, resistance (p.u.)
LN_X			= 5;	% x, reactance (p.u.)
LN_B			= 6;	% b, total line charging susceptance (p.u.)
LN_PF			= 7;	% real power injected at "from" bus end (MW)  
LN_QF			= 8;	% reactive power injected at "from" bus end (MVAr) 
LN_PT			= 9;	% real power injected at "to" bus end (MW)		
LN_QT			= 10;	% reactive power injected at "to" bus end (MVAr)	 
LN_FVM		    = 11;   % voltage magnitude at "from" bus 
LN_TVM		    = 12;   % voltage magnitude at "to" bus 
LN_OPTPF		= 13;	
LN_OPTQF		= 14;	
LN_OPTPT		= 15;		
LN_OPTQT		= 16;	
LN_OPTFVM	    = 17;
LN_OPTTVM	    = 18;
return;


