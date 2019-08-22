%% This demo recreates the numerical experiments with the random A
%
% REQUIREMENTS:
% Install UNLocBox from: https://epfl-lts2.github.io/unlocbox-html/
% Install SuiteLasso from: http://www.math.nus.edu.sg/~mattohkc/SuiteLasso.html


%creating the problem
[A,b,x,Lip] = createProblem(1000,0.02,1);
figure(1); plot(x); title('original'); ylim([-0.2,4.2]);


%setting parameters
it = 1000;
epsilon = 1e-4;
rho = 0;


%% solve problem using sparse Kaczmarz

lambda1 = 0.06804; %TV 
lambda2 = 0.055695; %l1

tic
[x_recSKA,con1,l] = sparseKaczmarz(A.',b,it,rho,lambda1, lambda2,0,0,epsilon);
tt1 = toc;
pp1 = psnr(x_recSKA,x);
ss1 = ssim(x_recSKA,x);
fprintf('sparse KA: PSNR=%6.4f, SSIM=%6.4f, TIME=%6.4f\n',pp1,ss1,tt1);

%plot convergence:
figure(2), semilogy(con1(1:l)); xlabel('Iterations'); ylabel('Relative Norm');

%% solve problem using SuiteLasso

%setting options:
opts.Lip=Lip;
opts.stoptol = epsilon;
opts.printyes = 0;

lambda1 = 0.003135;
lambda2 = 0.119966;

tic
x_recSL = Fused_Lasso_SSNAL(A,b,size(A,2),lambda1,lambda2,opts);      
tt2=toc;
pp2 = psnr(x_recSL,x);
ss2 = ssim(x_recSL,x);
fprintf('S Lasso  : PSNR=%6.4f, SSIM=%6.4f, TIME=%6.4f\n',pp2,ss2,tt2);

%% plotting the results:

figure(3), plot(x), hold on, plot(x_recSKA), plot(x_recSL); hold off;
legend('original','sparse KA','SSNAL');


%% greedy approach to find the best lambda1 and lambda2

    %initial values:
    lambda1 = 0.06;
    lambda2 = 0.05;
    ffun = @(x,y) sparseKaczmarz(A.',b,it,rho,x,y,0,0,epsilon);
    pval2 = psnr(ffun(lambda1,lambda2),x);
    [lambda1_sKA,lambda2_sKA,psnr_sKA] = refineLambda(lambda1,lambda2,ffun,pval2,x);

    %%
    %initial values
    lambda1 = 0.005;
    lambda2 = 0.1;
    ffun = @(x,y) Fused_Lasso_SSNAL(A,b,size(A,2),x,y,opts);   
    pval2 = psnr(ffun(lambda1,lambda2),x);
    [lambda1_SL,lambda2_SL,psnr_SL] = refineLambda(lambda1,lambda2,ffun,pval2,x);


