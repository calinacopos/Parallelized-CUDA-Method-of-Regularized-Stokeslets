% file to load
filename = 'pressure_shear.txt';

% extract the data
Z = load(filename);
x=reshape(Z(:,1),32,32);
y=reshape(Z(:,2),32,32);
z=reshape(Z(:,3),32,32);

% create 3D mesh grids
% mesh(x,y,z)
pcolor(x,y,z)
subplot(1,1,1);
axis square
shading flat
colorbar