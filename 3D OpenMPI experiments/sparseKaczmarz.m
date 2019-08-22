function [ x,con,l ] = sparseKaczmarz( A,b,iterations,rho,lambda1,lambda2,shuff,enforceReal,enforcePositive,epsilon,dim,x0 )
% sparse Kaczmarz algorithm for 3D
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
%             dim: Dimension of 3D x
%              x0: Starting value for iteration (optional)
%
% OUTPUT:
%               x: Reconstructed signal
%             con: Convergence
%               l: Number of iterations
%
% AUTHOR:
% Florian Lieb, August 2019


% initialization
[N, M] = size(A);
if nargin < 12
x = complex(zeros(N,1));
else 
    x=x0;
end
residual = complex(zeros(M,1));
rowIndexCycle = 1:M;
con = zeros(1,iterations);
x_k = x;

% calculate the energy of each frequency component
energy = vecnorm(A);

% may use a randomized Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
end

%thr = 1e-4;
st = @(s,T) s.*max(1 - T./abs(s),0);
ew = @(s,T) s.*max(1 - (T./abs(s)).^2,0);


% estimate regularization parameter
lambdZero = sum(energy.^2)/N;
lambdIter = rho*lambdZero; 

% settings for prox_tv3d
params = struct;
params.verbose = 0;
params.useGPU=0;

for l = 1:iterations
    
    for m = 1:M
        k = rowIndexCycle(m);
        
        if energy(k) > 0
            tmp = A(:,k).'*x;
            beta = (b(k) - tmp - sqrt(lambdIter)*residual(k)) / (energy(k)^2 + lambdIter);
            x = x + beta*conj(A(:,k));
            residual(k) = residual(k) + beta*sqrt(lambdIter);
        end
    end
    
    if enforceReal && ~isreal(x)
        x = complex(real(x),0);
    end
    
    if enforcePositive
        x(real(x) < 0) = 0;
    end
    
   
    x = prox_tv3d(reshape(x,dim'),lambda1,params);
    x = x(:);
    x = ew(x,lambda2);
    
    con(l) = norm(x_k - x,'fro')./(norm(x_k,'fro')+1e-3);
    
    if con(l) < epsilon
        break
    end
    
    x_k  = x;
end

end