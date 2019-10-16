function [Z_hat, E_hat, Q_hat, iter, L ] = NSLRG(X, lambda1, lambda2, lambda3, epsilon1, epsilon2, maxIter)
%% %%%%% L(Z,E,Q,Y1,Y2,mu) = |Z|_*  +lambda1 tr(ZLZ')+lambda2 |E|_1+ lambda3 * |Q|_1+ <Y1,X-XZ-E> +<Y2,Z-Q>+ mu/2 *( |X-XZ-E|_F^2+|Z-Q|_F^2), s.t. Z=Z',Z>0
%% %%%%%%
%% %% This matlab code used the Linearized Alternating Direction Method with Adaptive Penalty (LADMAP)
%%                       to solve non-negative symmetric low-rank representation graph regularized (NSLRG) method.
% %%----------------------------------------------------
%% %%-------------------------------------------------
%       X - m x n matrix of observations data (required input)
%

%% %%%  main code

addpath PROPACK;

if nargin < 2
    lambda1 = 100; 
end

if nargin < 3
    lambda2 = 100; 
end   

if nargin < 4
    lambda3 = 100; 
end   

if nargin < 5
   epsilon1 = 1e-6;
end

if nargin < 6
    epsilon2 = 1e-2; 
end 

if nargin < 7
    maxIter = 15; 
end 

%Initialize the Laplacian matrix L
     
 fea = X';      
 options = [];     
 options.NeighborMode = 'KNN';     
 options.k = 5;     
 options.WeightMode = 'HeatKernel';
 W = constructW(fea,options);
 

[mFea,nSmp]=size(X);
DCol = full(sum(W,2));
D = spdiags(DCol,0,nSmp,nSmp);
D = full(D);
L = D - W; 

% initialization
[m, n] = size(X);
Z_hat = zeros( n, n);
Q_hat = zeros( n, n);
E_hat = zeros( m, n);
Y1_hat = zeros( m, n);
Y2_hat = zeros( n, n);
rho = 2.5; % this one can be tuned
mu = 1e-3; % this one can be tuned
mu_max = 1e6; % this one can be tuned


iter = 0;
total_svd = 0;
converged = false;
stopCriterion = 1;
sv = 10;

while ~converged       
    iter = iter + 1;

   %updating Z
   Z_k = Z_hat;
   H = lambda1*(Z_hat*L' + Z_hat*L) + mu*(Z_hat - Q_hat + Y2_hat/mu) + mu*X'*(X*Z_hat - X + E_hat - Y1_hat/mu);
   p = 2*lambda1*norm(L,2) + mu*(1 + (norm(X,2))^2);
   
    if choosvd(n, sv) == 1
         [U, S, V] = lansvd(((Z_hat - H/p)+(Z_hat - H/p)')/2, sv, 'L');
       
    else
        [U, S, V] = svd(((Z_hat - H/p)+(Z_hat - H/p)')/2, 'econ');
    end
    diagS = diag(S);
    svp = length(find(diagS > 1/p));
    if svp < sv
        sv = min(svp + 1, n);
    else
        sv = min(svp + round(0.05*n), n);
    end
     Z_hat = U(:, 1:svp) * diag(diagS(1:svp) - 1/p) * V(:, 1:svp)';
     %Non-negative constraints      
     Z_hat(Z_hat<0)=0; 

%updating E
    E_k = E_hat;
    temp_E = X - X*Z_hat + Y1_hat/mu;
    E_hat = max(temp_E - lambda2/mu, 0) + min(temp_E + lambda2/mu, 0);  
%updating Q     
    Q_k = Q_hat;
    temp_Q = Z_hat + Y2_hat/mu;
    Q_hat = max(temp_Q - lambda3/mu, 0) + min(temp_Q + lambda3/mu, 0);
    Q_hat = max(Q_hat,0);
    
    total_svd = total_svd + 1;
    
 %updating Y1 Y2   
  Y1_hat = Y1_hat + mu*(X - X*Z_hat - E_hat);
  Y2_hat = Y2_hat + mu*(Z_hat - Q_hat);
  mu = min(mu_max,rho*mu);

     %% stop Criterion  

  b = norm(Z_hat - Z_k,2);
  c = norm(Q_hat - Q_k,2);
  d = norm(E_hat - E_k,2);
  a = [p*b, mu*c, mu*d];
  a = max(a);
   if  a<=epsilon2
        rho_hat = rho;
    else rho_hat = 1;
    end
  
    stopCriterion1 = norm(X - X*Z_hat - E_hat,2)/norm(X,2);
    e = [p*norm(Z_hat - Z_k,2),mu*norm(Q_hat - Q_k,2),mu*norm(E_hat - E_k,2)];
        stopCriterion2 = max(e);
        if stopCriterion1<epsilon1 && stopCriterion2<epsilon2
            converged = 1;
        end
        
         if mod( total_svd, 10) == 0
        disp(['#svd ' num2str(total_svd) ' r(Z) ' num2str(rank(Z_hat))...
            ' |E|_0 ' num2str(length(find(abs(E_hat)>0)))...
            ' stopCriterion ' num2str(stopCriterion)]);
         end    
        
         if ~converged && iter >=maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;       
         end
end
   
