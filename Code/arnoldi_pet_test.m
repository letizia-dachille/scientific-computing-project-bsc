function [x,err,mv,t] = arnoldi_pet_test(H,alpha,v,x,m,m1,beta,maxit,tol)
% function [x,err,mv,t] = arnoldi_pet_test(H,alpha,v,x,m,m1,beta,maxit,tol)
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  x: unit positive initial
%		  m: dimension of the Krylov subspace
%		  m1: frequency of using the extrapolation based on trace
%		  beta, maxit: parameters used to flip-flop between the PET 
%					   method and the Arnoldi type algorithm
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

% Set r = 1 and k = 1.
r = 1;
k = 1;

% Errors and matrix-vector products.
err = [];
mv = 0;

% Compute the number of dangling nodes l and the trace mu.
l = sum(dang);
mu = 1 + alpha * (l / n - 1);

tic;

while r > tol
   
    % Run 2 iterations of Arnoldi type algorithm.
	[x,r,mv] = short_ar(H,alpha,v,x,m,tol,mv);
	
    % If the residual norm satisfies the prescribed tolerance tol, then 
    % stop, else continue.
	if r <= tol
        break;
	end
    
    % Run the PET method with x as the initial guess, where x is the 
    % approximation vector obtained from the Arnoldi type algorithm.
    xn = x;
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
			
			if (mod(k,m1)==0)
                x = xn - (mu - 1) * x;
                x = x / norm(x,1);
                r = norm(x - xn,2);
                if r <= tol
                    break;
                end
                xn = x;
			end
        end
        x = xn;
        if r / r1 > beta 
           restart = restart + 1;
        end
    end
end

t = toc;