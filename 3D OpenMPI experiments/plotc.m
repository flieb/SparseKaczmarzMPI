function [] = plotc(c,dim,fig,filename)

x = reshape(c,dim');
mx = ceil(dim(1)/2);
my = ceil(dim(2)/2);
mz = ceil(dim(3)/2);

if nargin < 4
    filename = '';
end

ww = 800;
hh = 400;
hh2 = 420;



f=figure(fig);
%f.Position(1:2) = [680 280];

if ~strcmp(filename,'../concentrationPhantom3D.mdf')
    
    f.Position(3:4) = [ww hh];
    
    a = squeeze(x(:,:,mx));
    b = squeeze(x(:,my,:));
    c = squeeze(x(mz,:,:));
    maxval = max([a(:); b(:); c(:)]);
    minval = min([a(:); b(:); c(:)]);

    ax1 = subplot(131);
    imagesc(a,[minval maxval]);
    title(['z=' num2str(mz)]);
    axis image;
    axis off;
    ax1.XTick = []; ax1.YTick = [];
    ax1.Position(1) = 0.08;
    ax1.Position(3) = 0.3;

    ax2 = subplot(132);
    imagesc(b,[minval maxval]);
    title(['y=' num2str(my)]);
    axis image;
    axis off;
    ax2.XTick = []; ax2.YTick = [];
    daspect([1 0.5 0.5])
    ax2.Position(3) = ax1.Position(3)/2;

    ax3 = subplot(133);
    imagesc(c,[minval maxval]);
    title(['x=' num2str(mx)]);
    axis image;
    axis off;
    ax3.XTick = []; ax3.YTick = [];
    daspect([1 0.5 0.5])
    ax3.Position(3) = ax2.Position(3);
    ax3.Position(1) = ax2.Position(1)+ax2.Position(3)+(ax2.Position(1)-(ax1.Position(1)+ax1.Position(3)));
    
    cb = colorbar; 
    ylabel(cb,'mmol/l');
    ax3.Position(3)=ax2.Position(3);
else
    f.Position(3:4) = [ww hh2];
    
    ax = subplot(231);
    imagesc((x(:,:,7)));
    title('z=7');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];

    ax = subplot(232);
    imagesc(squeeze(x(:,7,:)));
    title('y=7');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];

    ax = subplot(233);
    imagesc(squeeze(x(7,:,:)));
    title('x=7');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];
    
    ax = subplot(234);
    imagesc((x(:,:,13)));
    title('z=13');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];

    ax = subplot(235);
    imagesc(squeeze(x(:,13,:)));
    title('y=13');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];

    ax = subplot(236);
    imagesc(squeeze(x(13,:,:)));
    title('x=13');
    axis image;
    axis off;
    ax.XTick = []; ax.YTick = [];
end