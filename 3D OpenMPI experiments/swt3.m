function T = swt3(X,level,wname)
% SWT3 Stationary Wavelet Transform in 3D
% Copyright (C) 2020 Florian Lieb

if level>2
    error('level>2 not implemented');
end
if level<1
    error('level must be 1 or 2');
end

[LoD,HiD,LoR,HiR] = wfilters(wname);

%init
dec = [];

%get 1st level undecimated DWT
T = myudwt3(X,LoD,HiD,LoR,HiR,'ppd');

if level==2
    %upsample filters
    LoD(2,length(LoD)) = 0;
    LoD = LoD(:)';
    HiD(2,length(HiD)) = 0;
    HiD = HiD(:)';
    
    %get undecimated DWT with upsampled filter
    TT = myudwt3(X,LoD,HiD,LoR,HiR,'ppd');
    
    %unwrap cell structure
    T.dec(:,:,:,9:16) = TT.dec;
end


end 

%-----------------------------------------------------------------------%
function wt = myudwt3(X,LoD,HiD,LoR,HiR,mode)
    sX = size(X);
    dec = zeros([sX 8]);

    LoD = LoD(:)';HiD=HiD(:)'; LoR=LoR(:)'; HiR=HiR(:)';
    permVect = [];
    
    %filter rows
    [a_Lo,d_Hi] = wdec1D(X,LoD,HiD,permVect,mode);
    
    %filter columns
    permVect = [2,1,3];
    [aa_Lo_Lo,da_Lo_Hi] = wdec1D(a_Lo,LoD,HiD,permVect,mode);
    [ad_Hi_Lo,dd_Hi_Hi] = wdec1D(d_Hi,LoD,HiD,permVect,mode);
    
    %filter 3rd dim
    permVect = [1,3,2];
    [dec(:,:,:,1),dec(:,:,:,2)] = wdec1D(aa_Lo_Lo,LoD,HiD,permVect,mode);
    [dec(:,:,:,5),dec(:,:,:,6)] = wdec1D(ad_Hi_Lo,LoD,HiD,permVect,mode);
    [dec(:,:,:,3),dec(:,:,:,4)] = wdec1D(da_Lo_Hi,LoD,HiD,permVect,mode);
    [dec(:,:,:,7),dec(:,:,:,8)] = wdec1D(dd_Hi_Hi,LoD,HiD,permVect,mode);
    
    %store results
    wt.sizeINI = sX;
    wt.filters.LoD = LoD;
    wt.filters.HiD = HiD;
    wt.filters.LoR = LoR;
    wt.filters.HiR = HiR;
    wt.mode = mode;
    wt.dec = dec;
end

%-----------------------------------------------------------------------%
function [L,H] = wdec1D(X,Lo,Hi,perm,dwtEXTM)

    if ~isempty(perm) , X = permute(X,perm); end
    sX = size(X);
    if length(sX)<3 , sX(3) = 1; end

    lf = length(Lo);
    lx = sX(2);
    lc = lx+lf-1;
    if lx<lf+1
        nbAdd = lf-lx+1;
        switch dwtEXTM
            case {'sym','symh','symw','asym','asymh','asymw','ppd'}
                Add = zeros(sX(1),nbAdd,sX(3));
                X = [Add , X , Add];
        end
    end

    switch dwtEXTM
        case 'zpd'             % Zero extension.

        %case {'sym','symh'}    % Symmetric extension (half-point).
        %    X = [X(:,lf-1:-1:1,:) , X , X(:,end:-1:end-lf+1,:)];

        case 'ppd'             % Periodized extension (1).
            X = [X(:,end-lf+2:end,:) , X , X(:,1:lf-1,:)];
    end

    L = convn(X,Lo,'same');
    H = convn(X,Hi,'same');

    switch dwtEXTM
        case 'zpd'
        case 'ppd'
            if length(Lo)==2
                L=L(:,1:sX(2),:);
                H=H(:,1:sX(2),:);
            else
                L=L(:,2:sX(2)+1,:);
                H=H(:,2:sX(2)+1,:);
            end
    end

    if ~isempty(perm)
        L = ipermute(L,perm);
        H = ipermute(H,perm);
    end
end


