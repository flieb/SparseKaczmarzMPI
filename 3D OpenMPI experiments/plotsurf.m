function [] = plotsurf(A,dim,fignum,factor,nlog)
%plots surface of SF in 3D mode
% Copyright (C) 2020 Florian Lieb

if nargin<5
    nlog = 0;
end

if (size(A)~=dim)
    A = reshape(A,dim');
end
if ~isreal(A)
    if nlog
        A = log10(abs(A));
    else
        A = abs(A);
    end
end


figure(fignum); cla, clf;
%subplot(131);

%factor = 0.8; %.22
ma = max(A(:));
mi = min(A(:));
dm = factor*(ma - mi);

%%
A = smooth3(A,'box',1);
A = permute(A,[3 2 1]);
patch(isocaps(A,factor),'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(A,factor),'FaceColor','blue','EdgeColor','none','FaceAlpha',1);
isonormals(A,p1);
%view(3);
view([-40 10]);
xlim([0 dim(1)]);
ylim([0 dim(2)]);
zlim([0 dim(3)]);
axis vis3d
ax = gca;
ax.ZDir = 'reverse';
ax.FontSize = 6;
%camlight left
camlight(-40,10);
%colormap('gray');
lighting gouraud
xlabel('z');
ylabel('y');
zlabel('x');
grid on
set(gcf, 'Color', 'w');
