function [ x,con,l,cngPS,cngSS,rel ] = regularizedKaczmarz( A,b,iterations,lambd,shuff,enforceReal,enforcePositive,epsilon,compfun,xrec,x0 )
% regularized Kaczmarz
% As published here: http://iopscience.iop.org/article/10.1088/0031-9155/55/6/003/meta on page 1582.
% Other references : Saunders M 1995: Solution of sparse rectangular
% systems using LSQR and CRAIG
% or
% From Dax A 1993: On row relaxation methods for large constrained 
% least squares problems
% Copyright (C) 2020 Florian Lieb

% initialization of the variable
[N, M] = size(A);
if nargin < 8
    epsilon = 1e-3;
end

if nargin > 8
    flag = true;
else
    flag = false;
end

if nargin < 11
x = complex(zeros(N,1));
else 
    x=x0;
end

if nargin <10
    xrec = zeros(size(x));
end

residual = complex(zeros(M,1));
rowIndexCycle = 1:M;

% calculate the energy of each frequency component
% energy = rowEnergy(A);
energy = sqrt(sum(abs(A.*A),1));

% may use a randomized Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
end

% estimate regularization parameter
lambdZero = sum(energy.^2)/N; %N
lambdIter = lambd*lambdZero; %lambdIter = lambd*lambdZero;
%lambdIter = 1e0;

con = zeros(1,iterations);
cngPS = con;
cngSS = con;
rel = con;
x_k = x;

for l = 1:iterations
    for m = 1:M
        k = rowIndexCycle(m);
        
        if energy(k) > 0
            tmp = A(:,k).'*x;
            beta = (b(k) - tmp - sqrt(lambdIter)*residual(k))./(energy(k).^2 + lambdIter);
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
    
    %inn = find(x==0);
    %A(inn,:) = 0;
    
    con(l) = norm(x_k - x,'fro')./(norm(x_k,'fro')+1e-3);
    if flag
        [cngPS(l),cngSS(l)] = compfun(x);
    end
    %cng(l) = norm(A.'*x-b,'fro')./norm(b,'fro');
    rel(l) = norm(x - xrec,'fro')./norm(xrec,'fro');
    x_k  = x;
    if con(l) < epsilon
        %fprintf('no. of it: %d\n',l);
        break
    end
end
end