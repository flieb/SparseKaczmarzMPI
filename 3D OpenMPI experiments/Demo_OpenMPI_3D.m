%% Script for OpenMPI Data
%
% REQUIREMENTS:
% Download 3D SystemMatrix and Phantom measurements from OpenMPI page:
% https://magneticparticleimaging.github.io/OpenMPIData.jl/latest/index.html
%
% in order to use the python toolbox PYWT, start Matlab from Anaconda
% prompt, needed for the SWT
% Otherwise choose other proximal operator in ska() or fistaMPI()


%%

clear, close all, clc;

filename = 'E:\_MPI2 Data\OpenMPI Data\3D\6.mdf';

dim = double(h5read(filename,'/calibration/size')); %should give [37 37 37]
snr = double(h5read(filename,'/calibration/snr'));
isBackgroundFrame = h5read(filename,'/measurement/isBackgroundFrame');
N = length(isBackgroundFrame);




%% load System Matrix, subtract BG and stack receive coils

%only load frequencies with SNR larger than 1.5:
snrval = 1.5;

S = [];
snr2 = [];

%loop over x, y, and z coil
for kk = 1:3
    
    indx = find(snr(:,kk)>snrval);
    snr2 = [snr2; snr(indx,kk)];

    % read the data (this might take a while)
    tmp = h5read(filename,'/measurement/data',[1 1 kk 1],[N size(snr,1) 1 1]);
    
    % convert to struct to complex
    tmp = complex(tmp.r,tmp.i);
    tmp = tmp(:,indx);
    
    M = size(tmp,2);
    
    % extract background frames
    BG2 = tmp(logical(isBackgroundFrame),:);
    S2 = tmp(~logical(isBackgroundFrame),:);
    clear tmp;
    
    % interpolate bg measurements
    BGe = interp1(double(0:dim(1):prod(dim)),BG2,double(1:prod(dim)),'pchip');
    
    % subtract background from measurements
    S = [S S2-BGe];
    clear S2 BGe
end

%% load phantom measurements

filename = '../shapePhantom3D.mdf';
%filename = '../resolutionPhantom3D.mdf';
%filename = '../concentrationPhantom3D.mdf';    


% load data
isBackgroundFrame = h5read(filename,'/measurement/isBackgroundFrame');
dataP = h5read(filename,'/measurement/data');
dataph = squeeze(double(dataP(:,:,:,~logical(isBackgroundFrame))));
dataph = mean(dataph(:,:,:),3);
databg = mean(squeeze(double(dataP(:,:,:,logical(isBackgroundFrame)))),3);
clear dataP

% convert to fourier domain
dd1 = fft(dataph,[],1);
up = dd1(1:26929,:);
dd1 = fft(databg,[],1);
ue = dd1(1:26929,:);
upe = [];
clear dd1 dataph isBackgroundFrame

% create frequency vector
freq_sf = linspace(0,1250000,26929)';
freq = [];

% restrict phantom measurements to snrval
for kk=1:3
    indx = find(snr(:,kk)>snrval);
    freq = [freq; freq_sf(indx)];
    upe = [upe; up(indx,kk)-ue(indx,kk)];
end
clear up ue

%% preconditioning the frequency rows

FreqVal = 70000; %Hz
SnrVal = 3;

ind_fourier = find(freq>=FreqVal & snr2>SnrVal);
A = S(:,ind_fourier);
b = upe(ind_fourier);


%% compute A^*A and its largest eigenvalue for FISTA 

% this might take quite some time, as AtA is very large (approx. 20GB)
tic
energy = 1./sqrt(sum(abs(A.*A),1));

A2=(A.*energy).';
b2 = b.*energy.'; 
Ap = A2';
AtA = A2'*A2;


no_it = 10;
[X,Y] = arnoldi(AtA,randn(prod(dim),1),no_it);
Z = Y(1:no_it,1:no_it);

% gamma should be around gamma = 0.004569
gamma = 1/eigs(Z,1); 
toc


%% our sparse Kaczmarz approach with NNG thresholding

    lambda = 2.5e-2;     %regularization parameter
    it = 1000;         %number of iterations
    random = 0;        %use randomized Kaczmarz
    epsilon = 1e-2;    %stopping criterion
    th_method = 'ew';  %thresholding, either 'ew' (NNG) or 'st' (soft-thresholding)

    tic
    c1_ew = ska(A,b,it,lambda,random,dim,epsilon,th_method);
    tt = toc;

    %correction factor 2^3 * 100 to get correct mmol/l: 
    corr_fac = 800;
    fignum = 1;
    surf_th = 0.001;

    %plot layer view
    plotc(c1_ew*corr_fac,dim,fignum); 
    fig = gcf; 
    fig.Name = 'SKA (NNG)';
    h = colorbar; 
    ylabel(h,'mmol/l');

    %plot 3D surf
    plotsurf(c1_ew,dim,fignum+1,surf_th);
    title(['SKA (NNG): ' num2str(tt)],'FontSize',10);
    daspect([ 0.5 1 0.5])


%% our sparse Kaczmarz approach with ST thresholding

    lambda = 3e-4;     %regularization parameter
    it = 1000;         %number of iterations
    random = 0;        %use randomized Kaczmarz
    epsilon = 1e-2;    %stopping criterion
    th_method = 'st';  %thresholding, either 'ew' (NNG) or 'st' (soft-thresholding)

    tic
    c1_st = ska(A,b,it,lambda,random,dim,epsilon,th_method);
    tt = toc;

    %correction factor 2^3 * 100 to get correct mmol/l: 
    corr_fac = 800;
    fignum = 3;
    surf_th = 0.001;

    %plot layer view
    plotc(c1_st*corr_fac,dim,fignum); 
    fig = gcf; 
    fig.Name = 'SKA (ST)';
    h = colorbar; 
    ylabel(h,'mmol/l');

    %plot 3D surf
    plotsurf(c1_st,dim,fignum+1,surf_th);
    title(['SKA (ST): ' num2str(tt)],'FontSize',10);
    daspect([ 0.5 1 0.5])
    

%% FISTA approach with NNG thresholding

    lambda = 5e-5;     %regularization parameter
    it = 1000;         %number of iterations
    epsilon = 1e-3;    %stopping criterion
    th_method = 'ew';  %thresholding, either 'ew' (NNG) or 'st' (soft-thresholding)

    tic
    c2_ew = fistaMPI(A2,b2,it,lambda,dim,epsilon,th_method,gamma,Ap);
    tt = toc;

    %correction factor 2^3 * 100 to get correct mmol/l: 
    corr_fac = 800;
    fignum = 5;
    surf_th = 0.001;
    
    %plot layer view
    plotc(c2_ew*corr_fac,dim,fignum); 
    fig = gcf; 
    fig.Name = 'FISTA (EW)';
    h = colorbar; 
    ylabel(h,'mmol/l');

    %plot 3D surf
    plotsurf(c2_ew,dim,fignum+1,surf_th);
    title(['Fista (EW) : ' num2str(tt)],'FontSize',10);
    daspect([ 0.5 1 0.5])
    
    
%% FISTA approach with NNG thresholding

    lambda = 1e-5;     %regularization parameter
    it = 1000;         %number of iterations
    epsilon = 1e-3;    %stopping criterion
    th_method = 'st';  %thresholding, either 'ew' (NNG) or 'st' (soft-thresholding)

    tic
    c2_ew = fistaMPI(A2,b2,it,lambda,dim,epsilon,th_method,gamma,Ap);
    tt = toc;

    %correction factor 2^3 * 100 to get correct mmol/l: 
    corr_fac = 800;
    fignum = 7;
    surf_th = 0.001;
    
    %plot layer view
    plotc(c2_ew*corr_fac,dim,fignum); 
    fig = gcf; 
    fig.Name = 'FISTA (ST)';
    h = colorbar; 
    ylabel(h,'mmol/l');

    %plot 3D surf
    plotsurf(c2_ew,dim,fignum+1,surf_th);
    title(['Fista (ST) : ' num2str(tt)],'FontSize',10);
    daspect([ 0.5 1 0.5])
    
%%  State-of-the-art Tikhonov approach, regularized Kaczmarz (reg. KA)

    lambda = 5e-3;     %regularization parameter
    it = 1000;         %number of iterations
    random = 0;        %use randomized Kaczmarz
    epsilon = 1e-2;    %stopping criterion

    tic
    c3 = regularizedKaczmarz(A,b,it,lambda,random,1,1,epsilon);
    tt = toc;

    %correction factor 2^3 * 100 to get correct mmol/l: 
    corr_fac = 800;
    fignum = 9;
    surf_th = 0.001;
    
    %plot layer view
    plotc(c3*800,dim,fignum); 
    fig = gcf; 
    fig.Name = 'Regularized Ka';
    h = colorbar; 
    ylabel(h,'mmol/l');

    %plot 3D surf
    plotsurf(c3,dim,fignum+1,surf_th);
    title(['Regularized Ka : ' num2str(tt)],'FontSize',10);
    daspect([ 0.5 1 0.5])