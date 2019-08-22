%% Script for OpenMPI Data
%
% REQUIREMENTS:
% Install UNLocBox from: https://epfl-lts2.github.io/unlocbox-html/
% Download 3D SystemMatrix and Phantom measurements from OpenMPI page:
% https://magneticparticleimaging.github.io/OpenMPIData.jl/latest/index.html




%% load System Matrix, subtract BG and stack receive coils


filename = '../SM_3D.mdf';
dim = double(h5read(filename,'/calibration/size'));
snr = double(h5read(filename,'/calibration/snr'));
isBackgroundFrame = h5read(filename,'/measurement/isBackgroundFrame');
framePermutation = h5read(filename,'/measurement/framePermutation');
isFramePermutation = h5read(filename,'/measurement/isFramePermutation');

data = h5read(filename,'/measurement/data');
data = double(complex(data.r,data.i));

BG = data(logical(isBackgroundFrame),:,:);
S = data(~logical(isBackgroundFrame),:,:);

clear data;

BGe = interp1(double(0:dim(1):prod(dim)),BG,double(1:prod(dim)),'pchip');
SBGC_ = S-BGe;

clear BG BGe S

SBGC_ = reshape(SBGC_,size(SBGC_,1),size(SBGC_,2)*3);

%% load phantom measurements

filename = '../shapePhantom3D.mdf';
%filename = 'resolutionPhantom3D.mdf';
%filename = 'concentrationPhantom3D.mdf';

isBackgroundFrame = h5read(filename,'/measurement/isBackgroundFrame');
dataP = h5read(filename,'/measurement/data');

dataph = mean(squeeze(double(dataP(:,:,:,~logical(isBackgroundFrame)))),3);
databg = mean(squeeze(double(dataP(:,:,:,logical(isBackgroundFrame)))),3);

clear dataP


dd1 = fft(dataph,[],1);
up = dd1(1:26929,:);

dd1 = fft(databg,[],1);
ue = dd1(1:26929,:);

clear dd1 dataph databg

%% preconditioning the frequency rows

freq_sf = linspace(0,1250000,26929)';
freq = repmat(freq_sf,3,1);
%bandpass filter: 70 - 1250kHz and SNR > 4:
ind_fourier = find(freq>=70000  & snr(:)>4);

%set final variables:
A = SBGC_(:,ind_fourier);

upe = up - ue;
b = upe(ind_fourier);

%% some parameter defintions

%Number of iterations:
it = 1000;
%Stopping criterion:
epsilon = 1e-3;


%% Regularized Kaczmarz

rho=50e-3;     

tic
[c0,con0,ll0] = regularizedKaczmarz(A,b,it,rho,0,1,1,epsilon);
tt = toc;


plotc((c0)*100,dim,1,filename); 
fig = gcf; 
fig.Name = 'Regularized Kaczmarz - Layer View';

plotsurf(c0,dim,2,0.05);
fig = gcf; 
fig.Name = 'Regularized Kaczmarz - Isosurface';
daspect([ 0.5 1 0.5])

fprintf('Regularized Kaczmarz : %5.3f s\n',tt);


%% sparse Kaczmarz

rho = 0e-15;
lambda1 = 0.06;
lambda2 = 0.01;

tic
[c1,con1,ll1]=sparseKaczmarz( A,b,it,rho,lambda1,lambda2,0,1,1,epsilon,dim);
tt = toc;


plotc((c1)*100,dim,3,filename); 
fig = gcf; 
fig.Name = 'Sparse Kaczmarz - Layer View';

plotsurf(c1,dim,4,0.05);
fig = gcf; 
fig.Name = 'Sparse Kaczmarz - Isosurface';
daspect([ 0.5 1 0.5])

fprintf('Sparse Kaczmarz : %5.3f s\n',tt);

%% convergence plot

figure(5), 
loglog(con0(1:ll0));
hold on,
loglog(con1(1:ll1));
hold off;
legend('reg. KA','sparse KA');
xlabel('Iterations');
ylabel('rel. norm');