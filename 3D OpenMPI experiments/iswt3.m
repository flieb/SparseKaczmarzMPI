function X = iswt3(W)
% ISWT3 Inverse Stationary Wavelet Transform in 3D
% Copyright (C) 2020 Florian Lieb


if size(W.dec,4)==8
    level = 1;
elseif size(W.dec,4)==16
    level=2;
else
    error('error');
end

[Nx,Ny,Nz] = size(W.dec(:,:,:,1));

LoR = W.filters.LoR;
HiR = W.filters.HiR;


Y4 = W.dec(:,:,:,end);

for jj=level:-1:1
    step = 2^(jj-1);
    last = step;
    for first = 1:last

        iRow = first:step:Nx;
        iCol = first:step:Ny;
        i3 = first:step:Nz;

        s = [length(iRow),length(iCol),length(i3)];

        sR = iRow(1:2:length(iRow));
        sC = iCol(1:2:length(iCol));
        s3 = i3(1:2:length(i3));

        perm = [1,3,2];
        Q = cell(4,1);
        for kk=1:4
            Q{kk} = wrec1D(W.dec(sR,sC,s3,2*(kk-1)+1),LoR,perm,s)...
                  + wrec1D(W.dec(sR,sC,s3,2*(kk-1)+2),HiR,perm,s); 
        end
%         Q{1} = wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 1),LoR,perm,s)...
%              + wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 2),HiR,perm,s);
%         Q{2} = wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 3),LoR,perm,s)...
%              + wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 4),HiR,perm,s); 
%         Q{3} = wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 5),LoR,perm,s)...
%              + wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 6),HiR,perm,s); 
%         Q{4} = wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 7),LoR,perm,s)...
%              + wrec1D(W.dec(sR,sC,s3,(jj-1)*8+ 8),HiR,perm,s); 


        perm = [2,1,3];
        Z={};
        Z{1} = wrec1D(Q{1},LoR,perm,s) + wrec1D(Q{2},HiR,perm,s);
        Z{2} = wrec1D(Q{3},LoR,perm,s) + wrec1D(Q{4},HiR,perm,s); 

        Y = wrec1D(Z{1},LoR,[],s) + wrec1D(Z{2},HiR,[],s);
        Y = circshift(Y,-1,3);
        Y = circshift(Y,-1,2);
        Y = circshift(Y,-1,1);


        sR = iRow(2:2:length(iRow));
        sC = iCol(2:2:length(iCol));
        s3 = i3(2:2:length(i3));

        perm = [1,3,2];
        Q = cell(4,1);
        for kk=1:4
            Q{kk} = wrec1D(W.dec(sR,sC,s3,2*(kk-1)+1),LoR,perm,s)...
                  + wrec1D(W.dec(sR,sC,s3,2*(kk-1)+2),HiR,perm,s); 
        end

        perm = [2,1,3];
        Z={};
        Z{1} = wrec1D(Q{1},LoR,perm,s) + wrec1D(Q{2},HiR,perm,s);
        Z{2} = wrec1D(Q{3},LoR,perm,s) + wrec1D(Q{4},HiR,perm,s); 


        Y2 = wrec1D(Z{1},LoR,[],s) + wrec1D(Z{2},HiR,[],s);

        Y4 = 0.5*(Y2 + Y);
        if jj==1
            W.dec(iRow,iCol,i3,1)=Y4;
        else
            W.dec(iRow,iCol,i3,9)= Y4;
        end
    
    end
end
%compnorm(Y4,X)
X = Y4;
end



function Y = wrec1D(X,R,perm,s)

    if ~isempty(perm)
        X = permute(X,perm);
        s = s(perm);
    end
    
    sX = size(X);
    Z = zeros(sX(1),2*sX(2)-1,sX(3));
    Z(:,1:2:end,:) = X;
    Y = convn(Z,R);

    sX = size(Y,2);
    F  = floor((sX-s)/2);
    C  = ceil((sX-s)/2);
    Y  = Y(:,1+F(2):end-C(2),:);

    if ~isempty(perm) , Y = permute(Y,perm); end
end