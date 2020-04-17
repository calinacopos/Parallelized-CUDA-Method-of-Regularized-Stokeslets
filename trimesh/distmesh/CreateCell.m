clear all; clf;

xL = 6;
yL = 6;
rc = 5;
dx = 0.5;%0.02985;%uniform-256 %0.059%uniform-64
ds = 0.45;

[xx,yy]=meshgrid(-xL:dx:xL,-yL:dx:yL);
dd=sqrt(xx.^2+yy.^2)-5.0; 
fd=@(p) dcircle(p,0,0,0.5);

% uniform mesh
[p,tri]=distmesh2d(@dmatrix,@huniform,ds,[-xL,-yL;xL,yL],[],xx,yy,dd);

% adaptive mesh
%fh=@(p) 0.01+0.4*(0.5-dcircle(p,0,0,0.1));
%[p,tri]=distmesh2d(fd,fh,ds,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);

% compute the edge connectivity information
%
X = p(:,1);
Y = p(:,2);
TR = TriRep(tri,X,Y);
E  = edges(TR);

% identify the boundary points -- note first and last point are the same
%
B = convhull(X,Y);
 
% form membrane at a distance of 0.1
%
XM = 1.01*X(B);
YM = 1.01*Y(B);

% plot everything beautifully :)
%
figure(1);
clf;
tr = trimesh(tri,X,Y,'color','black','linewidth',2.0);
set(tr,'color',[0.5 0.5 0.5]);
hold on; axis square;
for i=1:length(XM)
    plot([XM(i) X(B(i))],[YM(i) Y(B(i))],'linewidth',2.0);
end
scatter(X,Y,'k*','LineWidth',2.0);
mem = plot(XM,YM,'-k','linewidth',1.5);
set(mem,'color',[0.5 0.5 0.5]);
scatter(XM,YM,'rs','fillcolor','red','LineWidth',2.0);
set(gca, 'visible', 'off') ;

% print to file
%
fid = fopen('CellCirclePointsN447.txt','w');
format = '%.16f %.16f \n';
fprintf(fid,format,[X Y]');
fclose(fid);

fid = fopen('CellCircleMembraneN447.txt','w');
format = '%.16f %.16f \n';
fprintf(fid,format,[XM YM]');
fclose(fid);

fid = fopen('CellCircleTrianglesN447.txt','w');
format = '%u %u %u \n';
fprintf(fid,format,tri');
fclose(fid);

fid = fopen('CellCircleBoundaryN447.txt','w');
format = '%u\n';
fprintf(fid,format,B');
fclose(fid);