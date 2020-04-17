% DELAUNEY TRIANGULATION FOR UNIT CIRCLE
% Updated: 11/29/12
% Calina Copos


[xx,yy]=meshgrid(-1.1:0.1:1.1,-1.1:0.1:1.1);  				% Generate grid
dd=sqrt(xx.^2+yy.^2)-1;  		      				% d(x,y) at grid points
[p,t]=distmesh2d(@dmatrix,@huniform,0.2,[-1,-1;1,1],[],xx,yy,dd);
