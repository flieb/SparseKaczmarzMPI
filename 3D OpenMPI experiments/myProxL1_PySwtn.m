function [a] = myProxL1_PySwtn(y,lambda,dim,method)
% l1 proximaal operator, based on 3D swt in python library pywt

% Important:
% requires to load the python library pywt prior to usage
% py.importlib.import_module('pywt')


%setup parameters
level = uint8(2);
wname = 'haar';
gamma = 1;

%reshape input y to 3D array
y = reshape(y,dim');


%embedding of y in larger cube such that its size is next power of 2
nx = 2^nextpow2(dim(1));
ny = 2^nextpow2(dim(2));
nz = 2^nextpow2(dim(2));
Z = zeros(nx,ny,nz);
nxb = max(floor((nx-dim(1))/2),1);
nyb = max(floor((ny-dim(2))/2),1);
nzb = max(floor((nz-dim(3))/2),1);
Z(nxb:nxb+dim(1)-1,nyb:nyb+dim(2)-1,nzb:nzb+dim(3)-1) = y;

%define thresholding rule: st - softthresholding, ew - empirical wiener (NNG)
switch method
    case 'st'
        ffun = @(x) py.pywt.threshold(x,gamma*lambda,'soft');
    case 'ew'
        ffun = @(x) py.pywt.threshold(x,gamma*lambda,'garrote');
    case 'ht'
        ffun = @(x) py.pywt.threshold(x,gamma*lambda,'hard');
    otherwise %just needed for testing:
        ffun = @(x) 2*x;
end


P = py.pywt.swtn(Z,wname,level);

for kk=1:length(P) %is equal to the number of level
    
    %get list of python dict keys for indexing
    keys = py.list(P{kk}.keys());
    
    for jj=1:double(py.len(P{kk})) %loop over all dict elements (aaa,aad,ada,add,...)
        
        %extract corresponding 3D coefficients
        dd = P{kk}{keys{jj}};
        %apply thresholding
        P{kk}{keys{jj}} = ffun(dd)-dd; 
    end
end
%reconstruction using iswtn
a2 = Z + 1/gamma*double(py.pywt.iswtn(P,wname));


%take out the resulting values of cube
a2 = a2(nxb:nxb+dim(1)-1,nyb:nyb+dim(2)-1,nzb:nzb+dim(3)-1);
a = a2(:);

end

