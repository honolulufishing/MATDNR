function[Ybus,Yf,Yt] = makeYbus(bus,Lnbr,trsfm,shtc,shtr,pq,pv)
% return bus admittance matrix

[BUS_I,BUS_TYPE,PD,QD,GS,BS,VM,VA,BASE_KV,VMAX,VMIN,OPTVM,OPTVA] = idx_Bus;
[TRS_NUM,F_BUS,T_BUS,TRS_R,TRS_X,TAP,TAPMAX,TAPMIN,TAPSERIES,TAPSTATE,TRS_PF,TRS_QF,TRS_PT,TRS_QT,...
      TRS_FVM,TRS_TVM,OPTTAP,TRS_OPTPF,TRS_OPTQF,TRS_OPTPT,TRS_OPTQT,TRS_OPTFVM,TRS_OPTTVM] = idx_TRSFM;
[LN_NUM,F_BUS,T_BUS,LN_R,LN_X,LN_B,LN_PF,LN_QF,LN_PT,LN_QT,LN_FVM,LN_TVM,LN_OPTPF,LN_OPTQF,...
      LN_OPTPT,LN_OPTQT,LN_OPTFVM,LN_OPTTVM] = idx_LNBR;
[SHTC_BUS,SHTC_QC,QCMAX,QCMIN,SHTC_SERIES,SHTC_STATE,OPTQC] = idx_SHTC;
[SHTR_BUS,SHTR_QR,QRMAX,QRMIN,SHTR_SERIES,SHTR_STATE,OPTQR] = idx_SHTR;

% find the total numbers of all the buses and all the line branches 
nl = size(Lnbr, 1)+size(trsfm, 1);		% number of lines
nb = size(bus,1);                       % number of bus

YLn = [];Ytrs = [];Bc =[];bt = [];
if ~isempty(Lnbr)
    YLn = 1./ (Lnbr(:, LN_R) + j * Lnbr(:, LN_X)); % series admittance
    Bc = Lnbr(:, LN_B);							   % line charging susceptance
end
if ~isempty(trsfm)
    Ytrs = 1./ (trsfm(:, TRS_R) + j * trsfm(:, TRS_X));	% series admittance
    bt = trsfm(:,TAP);							    % consideration without  phase shifter
end

Ytt = [YLn + j*Bc/2; Ytrs];
Yff = [YLn + j*Bc/2; Ytrs./bt.^2];
Yft = - [YLn ; Ytrs./bt];
Ytf = - [YLn ; Ytrs./bt];

Y0=zeros(nb,1);  % Note that if we ignore the branches to the grounds, we have to use cap. to provide reactive power in power flow calculation.
%Y0 = bus(:,GS)+j*bus(:,BS);	                         % vector of shunt admittances

% consider the admittance of capacity and shunts  
%if size(shtc,1)>0
%   cbus = shtc(:,SHTC_BUS);
 %  Y0(cbus) = Y0(cbus) + j *shtc(:, SHTC_QC);	     % plus generation
%end
%if size(shtr,1)>0
%    rbus=shtr(:,SHTR_BUS);
%    Y0(rbus) =Y0(rbus) - j *shtr(:, SHTR_QR);	     % plus generation
%end



%% build Ybus
  if  ~isempty(trsfm)
    f = [Lnbr(:, F_BUS);trsfm(:,F_BUS)];			% list of "from" buses
    t = [Lnbr(:, T_BUS);trsfm(:,T_BUS)];			% list of "to" buses
  else
    f = Lnbr(:, F_BUS);
    t = Lnbr(:, T_BUS);
  end

	Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);		% connection matrix for line & from buses
	Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);		% connection matrix for line & to buses
	Ybus = spdiags(Y0, 0, nb, nb) + ...				% shunt admittance
        Cf * spdiags(Yff, 0, nl, nl) * Cf' + ...	% Yff term of branch admittance
		Cf * spdiags(Yft, 0, nl, nl) * Ct' + ...	% Yft term of branch admittance
		Ct * spdiags(Ytf, 0, nl, nl) * Cf' + ...	% Ytf term of branch admittance
		Ct * spdiags(Ytt, 0, nl, nl) * Ct';			% Ytt term of branch admittance
	
	i = [[1:nl]'; [1:nl]'];		                    % double set of row indices	
	Yf = sparse(i, [f; t], [Yff; Yft]);
	Yt = sparse(i, [f; t], [Ytf; Ytt]);
 
% Sparse Symmetric Matrix Figure
% clf;
% figure(1)
% spy(Ybus);
% title('Ybus Sparse Symmetric Matrix')
% nz = nnz(Ybus);
% pct = 100 / numel(Ybus);
% xlabel(sprintf('nonzeros=%d (%.3f%%)',nz,nz*pct));   
    
     
% % build Bpp
% Bpp = -imag(Ybus);
% % build Bp(考虑以下情况:不计节点对地电纳，不计支路对地电纳，不计变压器变比，不计支路电阻，不计变压器电阻)
% %%temp_bus(:, BS) = zeros(nb, 1);				%% zero out shunts at buses
% %%temp_Lnbr(:, LN_B) = zeros(nl, 1);		%% zero out line charging shunts
% %%temp_trsfm(:, TAP) = ones(nt, 1);			%% cancel out taps
% %%temp_Lnbr(:, LN_R) = zeros(nl, 1);		%% zero out line charging shunts
% %%temp_trsfm(:, TRS_R) = zeros(nt, 1);			%% cancel out taps
% YLn = 1./  (j * Lnbr(:, LN_X));	   
% Ytrs = 1./ ( j * trsfm(:, TRS_X));	%% series admittance
% tempYtt = [YLn ; Ytrs];%Yff = [YLn; Ytrs];
% tempYft = - [YLn ; Ytrs];%Ytf = - [YLn ; Ytrs];
% Y0 = bus(:,GS);
% tempYbus = spdiags(Y0, 0, nb, nb) + ...				%% shunt admittance
%     Cf * spdiags(tempYtt, 0, nl, nl) * Cf' + ...	%% Yff term of branch admittance
% 	Cf * spdiags(tempYft, 0, nl, nl) * Ct' + ...	%% Yft term of branch admittance
% 	Ct * spdiags(tempYft, 0, nl, nl) * Cf' + ...	%% Ytf term of branch admittance
% 	Ct * spdiags(tempYtt, 0, nl, nl) * Ct';			%% Ytt term of branch admittance
% Bp = -imag(tempYbus);
% 
% temp = Bp(:, [pv; pq])';
% newBp = temp(:, [pv; pq])';
% temp = Bpp(:, pq)';
% Bpp = temp(:, pq)';
% 
% %% factor B matrices
% [Lp, Up, Pp] = lu(newBp);
% [Lpp, Upp, Ppp] = lu(Bpp);
% Output Arguments: Lp,Up,Pp,Lpp,Upp,Ppp

