%e = sparse(e);
%v = sparse(v);

d = H * e;
dang = d==0;
dh = d + dang*n;
dh = 1./dh;

%dh = sparse(dh);
%A = alpha * (H + dang * e')' * diag(dh) + (1 - alpha) * v * e';

%y = A * x;
y = x .* dh;
y = H' * y + sum(dang .* y);
y = y * alpha + (1 - alpha) * v;

norm(x-y,2)