%% LOADING THE VARIABLES FOR THE TEST
load('data/generalized_brown.mat'); %different for each run
load('data/forcing_terms.mat');
% Problem dimension
% n= 2;
% n = 10^3;
n = 10^4;
% Starting point
% x0 = ones(n,1); % starting point banded trigonometric
% x0 = 1.2*ones(n,1); % starting point Rosenbrock 2d
% x0 = 10*ones(n,1); % different starting point
for i=1:n % starting point generalized brown
    if mod(i,2)==1
        x0(i)=0;
    else
        x0(i)=-1;
    end
end

alpha0 = 1;
rho = 0.8;
btmax = 50;
gamma = 0.9;
tolx = 1e-6;
tollgrad = 1e-7;
kmax = 10000;
c1= 1e-4;

% Select function and gradient to use
f = @generalized_brown;
gradf = @generalized_brown_grad;

%% RUN FLETCHER AND REEVES
disp('**** FLETCHER AND REEVES: START ****')
tic
[xk_fr, fk_fr, gradfk_norm_fr, k_fr, xseq_fr, btseq_fr] = FR_CG_bcktrck(x0, f, ...
    gradf, alpha0, kmax, tollgrad, c1, rho, btmax);
toc
disp('**** FLETCHER AND REEVES: FINISHED ******' )
disp('**** FLETCHER AND REEVES: RESULTS *****' )
disp('**************************************')
disp( ['f(xk_fr): ' , num2str(fk_fr)])
disp( ['grad(xk_fr): ' , num2str(norm(gradf(xk_fr),2))])
disp(['N. of Iterations: ', num2str(k_fr), '/', num2str(kmax), ';'])
disp('**************************************')

%% RUN INEXACT NEWTON METHOD
Hessf = @(x) 0;
FDgrad = 1; % exact gradient
FDHess = 'Jc'; %hessian approximated through the jacobian
pcg_maxit = 50;

fterms = fterms_suplin;
h =sqrt(eps)*norm(x0);
disp('**** INEXACT NEWTON METHOD: START ****')
tic
[xk_n, fk_n, gradfk_norm_n, k_n, xseq_n] = innewton_general(x0, f, gradf, Hessf, ...
    kmax, tollgrad, c1, rho, btmax, FDgrad, FDHess, h, fterms, pcg_maxit);
toc
disp('**** INEXACT NEWTON METHOD: FINISHED ******' )
disp('**** INEXACT NEWTON METHOD: RESULTS *****' )
disp('**************************************')
disp( ['f(xk_n): ' , num2str(fk_n)])
disp( ['grad(xk_n): ' , num2str(norm(gradf(xk_n),2))])
disp(['N. of Iterations: ', num2str(k_n), '/', num2str(kmax), ';'])
disp('**************************************')

%% RUN THE NELDER MEAD
disp('**** NELDER MEAD: START ****')
tic
[xk_nm, fk_nm, k_nm] = nelder_mead(x0, f, tolx);
toc
disp('**** NELDER MEAD: FINISHED ******' )
disp('**** NELDER MEAD: RESULTS *****' )
disp('**************************************')
disp( ['f(xk_nm): ' , num2str(fk_nm)])
disp( ['grad(xk_nm): ' , num2str(norm(gradf(xk_nm),2))])
disp(['N. of Iterations: ', num2str(k_nm), '/', num2str(kmax), ';'])
disp('**************************************')
