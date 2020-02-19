function [ x,con,l ] = ska( A,b,iterations,lambda,shuff,dim,epsilon,method)
%SKA  Sparse Kaczmarz based on undecimated wavelet tranform
%
%   Input parameters:
%       A          : System matrix
%       b          : Phantom measurement
%       iterations : Number of iterations
%       lambda     : Regularization parameter
%       shuff      : Flag, use randomized Kaczmarz
%       dim        : Dimension of x
%       epsilon    : Stopping criterion
%       method     : 'ew' or 'st' for NNG or soft-thresholding
%   Output parameters:
%       x          : Reconstructed tracer concentration
%       con        : Stopping criterion per iteration
%       l          : Number of iterations until stopping criterion
%
%
% Copyright (C) 2020 Florian Lieb

if nargin<8
    method = 'st';
end

if nargin <7
    method = 'st';
    epsilon = 1e-3;
end

%check dimensions
[N, M] = size(A);

if N~=prod(dim)
    error('dimension mismatch');
end

%set output:
x = complex(zeros(N,1));


%row energy:
energy = vecnorm(A);


% may use a randomized Kaczmarz
if shuff
    rowIndexCycle = randperm(M);
else
    rowIndexCycle = 1:M;
end


con = zeros(1,iterations);
x_k = x;

for l = 1:iterations
    for m = 1:M
        k = rowIndexCycle(m);
        x = x - (A(:,k).'*x - b(k))*conj(A(:,k))/energy(k)^2;
    end
    
    %choose the Proximal operator
    %x = myProxL1_swt3(real(x),lambda,dim,method); 
    x = myProxL1_PySwtn(real(x),lambda,dim,method); %this requires python

    %project to positive numbers
    x(x<0) = 0;
    
    %exclude zero rows from system matrix
    %inn = find(x==0);
    %A(inn,:) = 0;

    % stopping criterion
    con(l) = norm(x_k - x,'fro')./(norm(x_k,'fro')+1e-3); 
    x_k  = x;
    if con(l) < epsilon
        break
    end

end

x = real(x);
end
