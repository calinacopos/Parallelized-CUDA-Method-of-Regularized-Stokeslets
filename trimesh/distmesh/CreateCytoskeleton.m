xL = 2;
yL = 2;
rc = 1;
dx = 0.45; %0.2969; %0.059
ds = dx;

%[xx,yy]=meshgrid(-xL:dx:xL,-yL:dx:yL);
%dd=sqrt(xx.^2+yy.^2)-rc; 
%fd=@(p) dcircle(p,0,0,0.5);

%[p,tri]=distmesh2d(@dmatrix,@huniform,ds,[-xL,-yL;xL,yL],[],xx,yy,dd);
fd = @(p) sqrt(sum(p.^2,2))-9.999;
[p,tri]=distmesh2d(fd,@huniform,ds,[-9.999,-9.999;9.999,9.999],[]);

% compute the edge connectivity information
%
X = p(:,1); size(X)
Y = p(:,2);
TR = TriRep(tri,X,Y);
E  = edges(TR);

% identify the boundary points -- note first and last point are the same
%
B = convhull(X,Y);


%Xshift = X+15;
%Yshift = Y+15;
Xshift = X;
Yshift = Y;


fid = fopen('UnitCirclePointsN1788.txt','w');
format = '%e %e \n';
fprintf(fid,format,[Xshift Yshift]');
fclose(fid);

fid = fopen('UnitCircleTrianglesN1788.txt','w');
format = '%u %u %u \n';
fprintf(fid,format,tri');
fclose(fid);

fid = fopen('UnitCircleMedgEdgesN1788.txt','w');
format = '%u %u \n';
fprintf(fid,format,E');
fclose(fid);

fid = fopen('UnitCircleBoundaryN1788.txt','w');
format = '%u\n';
fprintf(fid,format,B');
fclose(fid);