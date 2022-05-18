% Input parameters.
% matrix: [1] wb-cs-stanford
%		  [2] web-Stanford
%		  [3] Stanford_Berkeley
%		  [4] web-Google
%		  [otherwise] random
% script: [1] power
%		  [2] pet
%		  [3] arnoldi_type_algorithm
%		  [4] thick_restarted_arnoldi
%		  [5] power_arnoldi
%		  [6] arnoldi_pet
%		  [7] arnoldi_pet_test
% alpha: damping factor
matrix = 1;
script = 6;
alpha = 0.997;

% Fixed parameters.
% beta: parameter used to flip-flop between the PET method and the thick 
%		restarted Arnoldi algorithm
% tol: prescribed tolerance
tol = 1e-8;
beta = alpha - 0.1;

% H: adjacency matrix
% m: dimension of the Krylov subspace
% p: number of approximate eigenpairs which are wanted
% maxit: parameter used to flip-flop between the PET method and the thick 
%		 restarted Arnoldi algorithm
% m1: frequency of using the extrapolation based on trace
switch matrix
	case 1
		H = load("../Test_matrices/wb-cs-stanford.mat").Problem.A;
		m = 5;
		p = 3;
		maxit = 8;
		m1 = 50;
		textm = 'wb-cs-stanford';
	case 2
		H = load("../Test_matrices/web-Stanford.mat").Problem.A;
		m = 5;
		p = 3;
		maxit = 12;
		m1 = 40;
		textm = 'web-Stanford';
	case 3
		H = load("../Test_matrices/Stanford_Berkeley.mat").Problem.A;
		m = 8;
		p = 5;
		maxit = 6;
		m1 = 40;
		textm = 'Stanford_Berkeley';
	case 4
		H = load("../Test_matrices/web-Google.mat").Problem.A;
		m = 9;
		p = 4;
		maxit = 10;
		m1 = 40;
		textm = 'web-Google';
	otherwise
		n = 100000;
		den = 1 / n;
		H = random_matrix(n,den);
		
		% -- Choose wisely --
		m = 5;
		p = 3;
		maxit = 8;
		m1 = 50;
		textm = 'random';
end

% Size of matrix.
n = size (H,1);

% Personalization vector.
e = ones (n,1);
v = e / n;

% x0: unit positive initial (norm 1)
% v1: unit positive initial (norm 2)
switch script
	case 1
		x0 = v;
		[x,err,mv,t] = power_method(H,alpha,v,x0,tol);
		texts = 'P';
	case 2
		x0 = v;
		[x,err,mv,t] = pet_method(H,alpha,v,x0,m1,tol);
		texts = 'PET';
	case 3
		v1 = v/norm(v,2);
		[x,err,mv,t] = arnoldi_type_algorithm(H,alpha,v,v1,m,tol);
		texts = 'AR';
	case 4
		v1 = v/norm(v,2);
		[x,err,mv,t] = thick_restarted_arnoldi(H,alpha,v,v1,m,p,tol);
		texts = 'TRA';
	case 5
		x0 = v;
		[x,err,mv,t] = power_arnoldi(H,alpha,v,x0,m,beta,maxit,p,tol);
		texts = 'P-A';
	case 6
		x0 = v;
		[x,err,mv,t] = arnoldi_pet(H,alpha,v,x0,m,m1,beta,maxit,p,tol);
		texts = 'A-PET';
	case 7
		x0 = v;
		[x,err,mv,t] = arnoldi_pet_test(H,alpha,v,x0,m,m1,beta,maxit,tol);
		texts = 'A-PET-T';
	case 8
		x0 = v;
		[x,err,mv,t] = power_arnoldi(H,alpha,v,x0,m,beta,maxit,p,tol);
		texts = 'P-A';
end

% Check.
check;

% Number of iterations
it = size(err,2);

% Print data on file
name = '../Txt/' + string(textm) + '_' + string(alpha) + '.txt';
file = fopen(name, 'a');
fprintf(file, "\nMatrix: %s\nMethod: %s\nAlpha: %s\n\n", textm, texts, string(alpha));
fprintf(file, "Number of iterations: %d\n", it);
fprintf(file, "Number of matrix-vector products: %d\n", mv);
fprintf(file, "Time elapsed: %s\n\n", string(t));
fclose(file);

error_plot(err, texts, alpha);
hold on