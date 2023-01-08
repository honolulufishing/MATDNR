
function [A_total,b_total, Aeq ,beq] = polyhedral_approx(x,id_x,n_k)
% Polyhedral Approximation of Conic Constraints
% n_k = 3;
% X = [x1,x2,x3, sigma, eta]
% x3^2 >= x1^2 + x2^2;

x1 = id_x(1,1);  
x2 = id_x(2,1);
x3 = id_x(3,1);

n_x = size(x,1);
n_var = 2*(n_k + 1) + n_x ;  % total variables

%---------------equality constraints------------%
% sigma_j = cos(pi/2^(j+1))sigma_(j-1) + sin(pi/2^(j+1))eta_(j-1)
A_eq = zeros(n_k,n_var);
for i = 1: n_k
    A_eq(i,n_x+1+i) = 1;
    A_eq(i,n_x+i) = -cos(pi/2^(i-1+2));
    A_eq(i,n_x+(n_k+1)+i) = -sin(pi/2^(i-1+2));
end
%----------------- END ----------------%

%---------------inequality constraints------------%
% -sigma^0 <=x1<= sigma^0
A_1 = zeros(2,n_var);
A_1(1,x1) = 1; A_1(1,n_x+1) = -1;
A_1(2,x1) = -1; A_1(2,n_x+1) = -1;



% -eta^0 <=x2<= eta^0
A_2 = zeros(2,n_var);
A_2(1,x2) = 1; A_2(1,n_x+(n_k+1)+1) = -1;
A_2(2,x2) = -1; A_2(2,n_x+(n_k+1)+1) = -1;

% eta^j >= -sin(pi/2^(j+1))sigma_(j-1) + cos(pi/2^(j+1))eta_(j-1)
A_3 = zeros(n_k,n_var);
for i = 1: n_k
    A_3(i,n_x+(n_k+1)+1+i) = -1;
    A_3(i,n_x+i) = -sin(pi/2^(i-1+2));
    A_3(i,n_x+(n_k+1)+i) = cos(pi/2^(i-1+2));
end

% eta^j >= sin(pi/2^(j+1))sigma_(j-1) - cos(pi/2^(j+1))eta_(j-1)
A_4 = zeros(n_k,n_var);
for i = 1: n_k
    A_4(i,n_x+(n_k+1)+1+i) = -1;
    A_4(i,n_x+i) = sin(pi/2^(i-1+2));
    A_4(i,n_x+(n_k+1)+i) = -cos(pi/2^(i-1+2));
end

% eta^n_k<=tan(pi/2^(n_k+1))*sigma^n_k
A_5 = zeros(1,n_var);
A_5(1,n_x+(n_k+1)+n_k+1) = 1;
A_5(1,n_x+(n_k+1)) = -tan(pi/2^((n_k+1)));

% sigma^n_k<=x3
A_6 = zeros(1,n_var);
A_6(1,n_x + n_k+1) = 1;
A_6(1,x3) = -1;
%----------------- END ----------------%

A_total = [A_1;
    A_2;
    A_3;
    A_4;
    A_5;
    A_6;
    ];

b_total = zeros(4 +2*n_k+2,1);

Aeq = A_eq;
beq = zeros(n_k,1);





