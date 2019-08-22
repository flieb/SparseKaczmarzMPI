function [ x,con,l,obb] = sparseKaczmarz( A,b,iterations,rho,lambda1,...
                                          lambda2,shuff,enforcePositive,...
                                          epsilon )
% sparse Kaczmarz algorithm
%
% INPUT:
%               A: Matrix
%               b: Measurement vector
%      iterations: Number of iterations
%             rho: Regularization parameter for ||x||_2
%         lambda1: Regularization parameter for ||x||_TV
%         lambda2: Regularization parameter for ||x||_1
%           shuff: Flag, if randomized Kaczmarz 
% enforcePositive: Flag, if projection onto positive numbers
%         epsilon: Stopping criterion
%
% OUTPUT:
%               x: Reconstructed signal
%             con: Convergence
%               l: Number of iterations
%             obb: Objective function
%
% AUTHOR:
% Florian Lieb, August 2019



% initializations
[N, M] = size(A);
x = complex(zeros(N,1));
residual = zeros(M,1);
rowIndexCycle = 1:M;
con = zeros(1,iterations);
obb = con;
x_k = x;

% calculate the rowenergy of A
energy = vecnorm(A);

% may use a randomized Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
end

%define some functionals:
st = @(s,T) s.*max(1 - T./abs(s),0);      % soft-thresholding
ew = @(s,T) s.*max(1 - (T./abs(s)).^2,0); % non-negative garrote
% objective functional
obj = @(x) 0.5*norm(A.'*x-b).^2 + rho*norm(x).^2 + ...
           lambda1*norm_tv1d(x) + lambda2*norm(x,1);


% estimate regularization parameter
lambdZero = sum(energy.^2)/N;
lambdIter = rho*lambdZero; 

% params for the tv-prox operator:
params = struct;
params.verbose = 0;
params.use_fast = 0; % set to 1 for Condats algorithm


for l = 1:iterations
    for m = 1:M
        k = rowIndexCycle(m);
        
        if energy(k) > 0
            tmp = A(:,k).'*x;
            beta = (b(k) - tmp - sqrt(lambdIter)*residual(k)) / ...
                   (energy(k)^2 + lambdIter);
            x = x + beta*A(:,k);
            residual(k) = residual(k) + beta*sqrt(lambdIter);
        end
    end
    
    
    if enforcePositive
        x(x < 0) = 0;
    end
    
    [x,info] = prox_tv1d(x,lambda1,params);
    x = ew(x,lambda2);

    %estimate convergence and objective functional    
    con(l) = norm(x_k - x)./(norm(x_k)+1e-3);
    %obb(l) = obj(x);
    
    %updates:
    x_k  = x;
    
    %check if converged
    if con(l) < epsilon
        break;
    end
end
end