function [a] = myProxL1_swt3(y,lambda,dim,method)
%% l1 proximaal operator, based on my own 3D swt implementation

%setup parameters
level = 2;
wname = 'haar';
gamma = 1;

%reshape input y to 3D array
y = reshape(y,dim');

ew  = @(s,T) s.*max(1-(T./abs(s)).^2,0);

%embedding of y in larger cube such that its size is next power of 2
nx = 2^nextpow2(dim(1));
ny = 2^nextpow2(dim(2));
nz = 2^nextpow2(dim(2));
nx=40;ny=nx;nz=40;
Z = zeros(nx,ny,nz);
nxb = max(floor((nx-dim(1))/2),1);
nyb = max(floor((ny-dim(2))/2),1);
nzb = max(floor((nz-dim(3))/2),1);
Z(nxb:nxb+dim(1)-1,nyb:nyb+dim(2)-1,nzb:nzb+dim(3)-1) = y;

%define thresholding rule: st - softthresholding, ew - empirical wiener (NNG)
switch method
    case 'st'
        ffun = @(x) wthresh(x,'s',gamma*lambda);
    case 'ew'
        ffun = @(x) ew(x,gamma*lambda);
    case 'ht'
        ffun = @(x) wthresh(x,'h',gamma*lambda);
    otherwise %just needed for testing:
        ffun = @(x) 2*x;
end

%swtn requires pywavelet package (py.importlib.import_module('pywt'))
P = swt3(Z,level,'haar');

P.dec = ffun(P.dec)-P.dec;


%reconstruction using iswtn
a2 = Z + 1/gamma*double(iswt3(P));


%take out the resulting values of cube
a2 = a2(nxb:nxb+dim(1)-1,nyb:nyb+dim(2)-1,nzb:nzb+dim(3)-1);
a = a2(:);

end

