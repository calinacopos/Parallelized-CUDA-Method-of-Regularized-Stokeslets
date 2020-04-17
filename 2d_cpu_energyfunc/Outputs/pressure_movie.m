clf; clear all;

% file to load
filename = 'pressure.txt';
imfile = 'testimage.tiff';

% extract the data
Z = load(filename);
N = 1;

x = reshape(Z(:,1),640,640);
y = reshape(Z(:,2),640,640);
p = reshape(Z(:,3),640,640);

for j=1:N 
    
%x=Z(1:256,1);
%y=Z(1:256,2);
%z=Z(1:256,3);

% create 3D mesh grids
% mesh(x,y,z)
    pcolor(xplot,yplot,pplot);
    axis square;
    shading interp;
    colorbar;
    M(j) = getframe;
%subplot(1,1,1);
end
movie(M)


% if making a moving, dump the plot to an image file
%
%saveas(gcf,imfile,'tiff');
