function [x,err,mv,t] = power_method(H,alpha,v,x,tol)
% function [x,err,mv,t] = power_method(H,alpha,v,x,tol)
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  x: unit positive initial
%		  tol: prescribed tolerance
% Output - x: PageRank vector
%		   err: vector containing errors
%		   mv: number of matrix-vector products
%		   t: time elapsed

% Size of matrix.
n = size (H,1);

% Matrix D.
e = ones (n,1);
d = H * e;
dang = d==0;
dh = d + dang*n;
dh = 1./dh;

% Set r = 1 and k = 0.
r = 1;
k = 0;

% Errors and matrix-vector products.
err = [];
mv = 0;

tic;

while r > tol

	% xn = A * x;
	xn = x .* dh;
	xn = H' * xn + sum(dang .* xn);
	xn = xn * alpha + (1 - alpha) * v;
	mv = mv + 1;

	k = k + 1;

	r = norm(xn - x,2);
	err = [err r];%#ok<AGROW>
	xn = xn / norm(xn,1);
	if r <= tol
		break;
	end
	
	x = xn;
end

t = toc;