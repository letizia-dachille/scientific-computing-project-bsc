function [x,err,mv,t] = power_arnoldi(H,alpha,v,x,m,beta,maxit,p,tol)
% function [x,err,mv,t] = power_arnoldi(H,alpha,v,x,m,beta,maxit,p,tol);
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  x: unit positive initial
%		  m: dimension of the Krylov subspace
%		  beta, maxit: parameters used to flip-flop between the PET 
%					   method and the thick restarted Arnoldi algorithm
%		  p: number of approximate eigenpairs which are wanted
%		  tol: prescribed tolerance
% Output - x: PageRank vector
%		   err: vector containing errors
%		   mv: number of matrix-vector products
%		   t: time elapsed

% Size of matrix.
n = size (H,1);

% Dangling nodes.
e = ones (n,1);
d = H * e;
dang = d==0;
dh = d + dang*n;
dh = 1./dh;

% Set r = 1, r0 = r and k = 1.
r = 1;
r0 = r;
k = 1;

% Errors and matrix-vector products.
err = [];
mv = 0;

% First run of thick restarted Arnoldi.
fr = 1;

% Parameter used to flip-flop between the PET method and the thick 
% restarted Arnoldi algorithm.
restart = 0;

tic;

% Run the power method once.
xn = x;
while r > tol && restart < maxit
    x = xn;

    % xn = A * x;
    xn = x .* dh;
    xn = H' * xn + sum(dang .* xn);
    xn = xn * alpha + (1 - alpha) * v;
	mv = mv + 1;
	
	k = k + 1;
    r = norm(xn - x,2);
    xn = xn / norm(xn,1);
	err = [err r]; %#ok<AGROW>
    if r < tol
        break;
    end
    if r / r0 > beta
        restart = restart + 1;
    end
    r0 = r;
end

x = xn;

while r > tol
	
    % Run 2 iterations of thick restarted arnoldi. Iterate
	% steps 2-8 of Algorithm 4 for the first run and steps 3-8 otherwise.
	if fr == 1
		[x,x1t,r,mv] = short_tra_fr(H,alpha,v,x,m,p,tol,mv);
		fr = 0;
	else
		[x,x1t,r,mv] = short_tra(H,alpha,v,x,m,p,tol,mv);
	end
    
    % If the residual norm satisfies the prescribed tolerance tol, then 
    % stop, else continue.
    if r <= tol
        break;
    end
    
    % Run the power method with x1t as the initial guess, where x1t is the 
    % approximation vector obtained from the thick restarted Arnoldi 
    % algorithm.
	x1t = abs(x1t);
    xn = x1t / norm (x1t,1);
    restart = 0;
    
    while r > tol && restart < maxit
        ratio = 0;
        r0 = r;
        r1 = r;
        while r > tol && ratio < beta
            x = xn;
            
            % xn = A * x;
            xn = x .* dh;
            xn = H' * xn + sum(dang .* xn);
            xn = xn * alpha + (1 - alpha) * v;
			mv = mv + 1;
            
            r = norm(xn - x,2);
            xn = xn / norm(xn,1);
            ratio = r / r0;
            r0 = r;
            k = k + 1;
			err = [err r];%#ok<AGROW>
        end
        x = xn;
        if r / r1 > beta 
           restart = restart + 1;
        end
    end
end

t = toc;