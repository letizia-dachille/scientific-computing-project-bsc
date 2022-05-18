function [x,err,mv,t] = arnoldi_type_algorithm(H,alpha,v,v1,m,tol)
% function [x,err,mv,t] = arnoldi_type_algorithm(H,alpha,v,v1,m,tol)
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  v1: unit positive initial
%		  m: dimension of the Krylov subspace
%		  tol: prescribed tolerance
% Output - x: PageRank vector
%		   err: vector containing errors
%		   mv: number of matrix-vector products
%		   t: time elapsed

% Errors and matrix-vector products.
err = [];
r = 1;
mv = 0;

tic;

while r > tol
	
    % Apply Algorithm 2 to form Vm1 = V_{m+1}, H_m, Hb = \bar{H_m}. 
    [Hb,Vm1] = arnoldi_process(H,alpha,v,v1,m);
	mv = mv + m;
    Vm = Vm1(:,1:m);
    
    Hb = Hb - [eye(m);zeros(1,m)];
    [~,~,V] = svd(Hb);
    v1 = Vm * V(:,m);
	
    r = svds(Hb,1,'smallestnz');
	err = [err r];%#ok<AGROW>
end

t = toc;

x = sign(sum(v1))*v1;
x = x/norm(x,1);