%rc = 9.0625;
% rc = 9;
% Aiming for 134 points on the boundary with 32 Lagrangian points
xL = 1.25;
yL = 1.25;
rc = 1;
dx = 2.95/32;
ds = dx*0.5;
% more refined
%[xx,yy]=meshgrid(-10:0.05:10,-10:0.05:10);
[xx,yy]=meshgrid(-xL:dx:xL,-yL:dx:yL);
dd=sqrt(xx.^2+yy.^2)-rc; 
fd=@(p) dcircle(p,0,0,0.5);

% define region of finer discretization
fh=@(p) 0.05+0.5*(0.5-dcircle(p,0,0,0.1));

% fd=inline('sqrt(sum((p-15).^2,2))-9.0625','p');
[p,tri]=distmesh2d(@dmatrix,@huniform,ds,[-xL,-yL;xL,yL],[],xx,yy,dd);
%[p,tri]=distmesh2d(fd,fh,0.02,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);

% compute the edge connectivity information
%
X = p(:,1);
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


fid = fopen('UnitCirclePointsN32.txt','w');
format = '%e %e \n';
fprintf(fid,format,[Xshift Yshift]');
fclose(fid);

fid = fopen('UnitCircleTrianglesN32.txt','w');
format = '%u %u %u \n';
fprintf(fid,format,tri');
fclose(fid);

fid = fopen('UnitCircleMedgEdgesN32.txt','w');
format = '%u %u \n';
fprintf(fid,format,E');
fclose(fid);

%BoundaryShift = 15+B;

fid = fopen('UnitCirlceBoundaryN32.txt','w');
format = '%e %e \n';
fprintf(fid,format,B');
fclose(fid);