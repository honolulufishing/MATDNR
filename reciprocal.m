function Y = reciprocal(X)
%RECIPROCAL  Invert the nonzero entries of a matrix elementwise.
%   Y = RECIPROCAL(X) has the same sparsity pattern as X
%	 (except possibly for underflow).

[m, n]  = size(X);
[i,j,Y] = find(X);
Y = sparse(i,j,1./Y,m,n);
Y = min(Y,1e8);
return;