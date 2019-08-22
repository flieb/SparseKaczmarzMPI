function [A,b,x,Lip] = createProblem(m,sigma,rownormalize)

% this function creates the problem
%
% INPUT:
%            m: number of rows
%        sigma: amount of noise
% rownormalize: normalize rows of A
%
% OUTPUT:
%            A: random matrix of size m x 2000
%            b: measurement vector
%            x: sparse signal of length 2000
%
% Author: 
% Florian Lieb, August 2019
%

if nargin<3
    rownormalize=1;
end

rng(0);

%fix n to 2000
n = 2000;


A = randn(m,n);

%sparse signal:
x = zeros(n,1);
x(200:350) = 2; 
x(750:850) = 1; 
x(1250:1500) = 4; 
x(1790:1795) = 2.8; 

%figure(1), plot(x)

%forward simulation:
tmp = A*x;

%adding noise:
noise = (sigma*std(tmp))^2*randn(size(tmp));
b = tmp + noise;

%computing Lipschitz constant for SSNAL algorithm
Amap = @(x) A*x;
ATmap = @(x) A'*x;
AATmap = @(x) Amap(ATmap(x));
eigsopt.issym = 1;
Lip = eigs(AATmap,length(b),1,'LA',eigsopt);
%figure(2); plot(tmp), hold on, plot(b,'--'), hold off;

if rownormalize
    energy = 1./vecnorm(A,2,2);
    A = A.*energy;
    b = b.*energy;
end