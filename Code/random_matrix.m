function A = random_matrix(n,d)
% function A = random_matrix(n,d)
% Input - n: size
%		  d: density
% Output - A: random, n-by-n, boolean matrix

A = sprand(n,n,d);
A = A ~= 0;