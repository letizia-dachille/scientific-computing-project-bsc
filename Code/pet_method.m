function [x,err,mv,t] = pet_method(H,alpha,v,x,m1,tol)
% function [x,err,mv,t] = pet_method(H,alpha,v,x,m1,tol)
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  x: unit positive initial
%		  m1: frequency of using the extrapolation based on trace
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

% Compute the number of dangling nodes l and the trace mu.
l = sum(dang);
mu = 1 + alpha * (l / n - 1);

% Errors and matrix-vector products.
err = [];
mv = 0;

tic;

xn = x;
while r > tol
    
    % Run the power iteration m_1 steps to obtain x_{m_1-1} and x_{m_1}.
    for i = 1:m1 
        x = xn;
        
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
    end
    
    if r <= tol
        x = xn;
        break;
    end
    
    % Use the extrapolation scheme based on x_{m_1-1}, x_{m_1} and mu.
    x = xn - (mu - 1) * x;
    x = x / norm(x,1);
    r = norm(x - xn,2);
    xn = x;
end

t = toc;