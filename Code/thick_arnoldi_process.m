function [Hb,V] = thick_arnoldi_process(H,alpha,v,Vp,Hp,m)
% function [Hb,V] = thick_arnoldi_process(H,alpha,v,Vp,Hp,m)
% Input - H: adjacency matrix
%		  alpha: damping factor
%		  v: personalization vector
%		  Vp: first part of the output V
%		  Hp: first part of the output Hb
%		  m: dimension of the Krylov subspace
% Output - Hb: upper Hessenberg matrix
%		   V: orthonormal basis of the Krylov subspace

% Size of matrix.
n = size (H,1);

% Calculate matrix D.
e = ones (n,1);
d = H * e;
dang = d==0;
dh = d + dang*n;
dh = 1./dh;

% Hb = \bar{H_m}, V = V_{m+1}.
Hb = zeros(m+1 , m);
Hb(1:size(Hp,1),1:size(Hp,2))= Hp;
V = zeros(n , m+1);
V(1:size(Vp,1),1:size(Vp,2))= Vp;

k = size(Vp,2);

for j = k : m
    
    % q = A * V(:,j);
    q = V(:,j) .* dh;
    q = H' * q + sum(dang .* q);
    q = q * alpha + (1 - alpha) * v * sum(V(:,j));
    
    for i = 1 : j
        Hb(i,j) = V(:,i)' * q;
        q = q - Hb(i,j) * V(:,i);
    end
    
    Hb(j+1,j) = norm(q,2);
    
    if Hb(j+1,j) == 0
        break;
    end
    
    V(:,j+1) = q / Hb(j+1,j);
end
