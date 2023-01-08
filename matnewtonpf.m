function [V,converged,k] = matnewtonpf(Ybus,V0,ref,pv,pq,Sbus,tol,pfmxt)
% Newton-Raphson Power flow Solver
% [V,Converged,K] = matnewtonpf(Ybus,V0,ref,pv,pq,Sbus,tol,pfmxt) returns the final voltage magnitude,
% phase angle and line flow and loss .
%
%  Computes partial derivatives of power injection with respect to voltage
%  magnitude and voltage angle in order to get a jacobian matrix.
%
%   S = diag(V) * conj(Ibus) = diag(conj(Ibus)) * V
%
%   Partials of V & Ibus w.r.t. voltage magnitudes
%       dV/dVm = diag(V./abs(V))
%       dI/dVm = Ybus * dV/dVm = Ybus * diag(V./abs(V))
%
%   Partials of V & Ibus w.r.t. voltage angles
%       dV/dVa = j * diag(V)
%       dI/dVa = Ybus * dV/dVa = Ybus * j * diag(V)
%
%   Partials of S w.r.t. voltage magnitudes
%       dS/dVm = diag(V) * conj(dI/dVm) + diag(conj(Ibus)) * dV/dVm
%              = diag(V) * conj(Ybus * diag(V./abs(V)))
%                                       + conj(diag(Ibus)) * diag(V./abs(V))
%
%   Partials of S w.r.t. voltage angles
%       dS/dVa = diag(V) * conj(dI/dVa) + diag(conj(Ibus)) * dV/dVa
%              = diag(V) * conj(Ybus * j * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = -j * diag(V) * conj(Ybus * diag(V))
%                                       + conj(diag(Ibus)) * j * diag(V)
%              = j * diag(V) * conj(diag(Ibus) - Ybus * diag(V))
%
% Note: In this program-package, we have considered all the shunts and capacity as constant loads.

% Author(s):Fang Liu, Chao Lei
% $Date:2012/05/04 17:30$


% Step 1: Initialize
j = sqrt(-1);
converged = 0;
V = V0;
Va = angle(V);
Vm = abs(V);

% Step 2: Set up indexing for updating V
npv	= length(pv);
npq	= length(pq);
n = npv + npq + 1;
j1 = 1;			     j2 = npv;			       % j1:j2 - V angle of pv buses
j3 = npv + 1;	     j4 = npv + npq;		   % j3:j4 - V angle of pq buses
j5 = npv + npq + 1;	 j6 = npv + npq + npq;     % j5:j6 - V mag of pq buses

% Step 3: Examine the initial values
% calculate power mismatches
mis = V .* conj(Ybus * V) - Sbus;
F = [	real(mis([pv; pq]));
    imag(mis(pq))	];

% check tolerance
normF = norm(F, inf);
if normF < tol
    converged = 1;
end

%??--------------------- Do iteration loop ---------------------?%
k = 1;
while  k < pfmxt

    % Step 4: Calculate power mismatches
    mis = V .* conj(Ybus * V) - Sbus;
    F = [	real(mis(pv));
        real(mis(pq));
        imag(mis(pq))	];

    % Step 5: Calculate jacobian elements
    Ibus = Ybus * V;
    if issparse(Ybus)			%% sparse version (if Ybus is sparse)
        diagV		= spdiags(V, 0, n, n);
        diagIbus	= spdiags(Ibus, 0, n, n);
        diagVnorm	= spdiags(V./abs(V), 0, n, n);
    else						%% dense version
        diagV		= diag(V);
        diagIbus	= diag(Ibus);
        diagVnorm	= diag(V./abs(V));
    end

    dSbus_dVm = diagV * conj(Ybus * diagVnorm) + conj(diagIbus) * diagVnorm;
    dSbus_dVa = j * diagV * conj(diagIbus - Ybus * diagV);

    % j11(angle) is a (n-1)¡Á(n-1) matrix
    % j12(voltage) is a (n-1)¡Á(nPQ) matrix
    % j21(angle) is a (nPQ)¡Á(n-1) matrix
    % j22(voltage) is a (nPQ)¡Á(nPQ) matrix
    j11 = real(dSbus_dVa([pv; pq], [pv; pq]));
    j12 = real(dSbus_dVm([pv; pq], pq));
    j21 = imag(dSbus_dVa(pq, [pv; pq]));
    j22 = imag(dSbus_dVm(pq, pq));

    % Jacobian matrix
    J = [	j11 j12;
        j21 j22;	];

    % Step 6: Compute update step
    dx = -(J \ F);

    % Step 7: Update voltage
    if size(pv,1)>0,
        Va(pv) = Va(pv) + dx(j1:j2);
    end
    Va(pq) = Va(pq) + dx(j3:j4);
    Vm(pq) = Vm(pq) + dx(j5:j6);
    V = Vm .* exp(j * Va);

    % Record the increments
    iterbus = zeros(n-1,4);
    iterbus(:,1) = dx(j1:j4); % PV + PQ: angle
    iterbus(1:npq,2) = dx(j5:j6);       % PV: magnitude
    iterbus(:,3) = F(1:npv+npq,1);      % real power
    iterbus(1:npq,4) = F(n:n+npq-1,1);  % reactive power

    [row,column] = size(iterbus);
    for i = 1:column
        iterbuscell{k,i} = iterbus(:,i);
    end
    
    % Step 8: Check for convergence
    mis = V .* conj(Ybus * V) - Sbus;
    F = [	real(mis(pv));
        real(mis(pq));
        imag(mis(pq))	];
    normF = norm(F, inf);
    if normF < tol
        converged = 1;
        break;
    end

    % k equals kmax or not
    if (k == pfmxt);
        converged = 0;
        disp('Number of iterations exceeded');
        break;
    end

    % Increase iterations
    k = k + 1;

end

% figure(2)
% % Show the graph of algorithm convergence
% realFrd = sum((cat(2,iterbuscell{:,3})),1);
% imagFrd = sum((cat(2,iterbuscell{:,4})),1);
% % Show the graph of algorithm convergence
% xk = [1:k];
% drealFrd = sqrt(realFrd.^2);
% dimagFrd = sqrt(imagFrd.^2);
% plot(xk,drealFrd,'k',xk,dimagFrd,'k','LineWidth',2)
% title('The Convergence Speed of Newton-Raphson Algorithm')
% xlabel('Iterations');
% ylabel('Speed');
% legend('Real Mismatches','Reactive Mismatches')
% grid on