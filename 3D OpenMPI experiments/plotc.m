function [maxval] = plotc(c,dim,fig)
% plots tracer concentration in 3D layer view

if nargin<5
    mv = 0;
end

x = reshape(c,dim');
mx = ceil(dim(1)/2);
my = ceil(dim(2)/2);
mz = ceil(dim(3)/2);

if nargin < 4
    filename = '';
end

ww = 1000;
hh = 365;
hh2 = 320;



f=figure(fig);
clf
%f.Position(1:2) = [680 280];

    
f.Position(3:4) = [hh ww];

a = squeeze(x(:,:,mx));
b = squeeze(x(:,my,:));
c = squeeze(x(mz,:,:));

maxval = max([a(:); b(:); c(:)]);
maxval = max(maxval,mv);

minval = min([a(:); b(:); c(:)]);

if maxval-minval==0
    maxval=1;
end

ax=subplot(311);
imagesc(fliplr(a.'),[minval maxval]);
axis image, axis off;
daspect([1 1 1])

ax=subplot(312);
imagesc(fliplr(b.'),[minval maxval]);
axis image, axis off;
daspect([1 1 1])

ax=subplot(313);
imagesc(c.',[minval maxval]);
axis image, axis off;
daspect([1 1 1])
