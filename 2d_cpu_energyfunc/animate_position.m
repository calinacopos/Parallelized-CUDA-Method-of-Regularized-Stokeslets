% 
% animation of position of confined prescribed blebbing
% 

clear all; clf;

load all.txt;
Npts = 64;
Nframes = length(all)/Npts;
Frames = moviein(Nframes);
speed = 10;

for i = 1:Nframes
    startPos = (i-1)*Npts + 1; 
    endPos = i*Npts;
    
    t = all(startPos:endPos,1);
    x = all(startPos:endPos,2);
    y = all(startPos:endPos,3);
    fx = all(startPos:endPos,4);
    fy = all(startPos:endPos,5);
    vx = all(startPos:endPos,6);
    vy = all(startPos:endPos,7);
     
    % Figure 1
    figure(1);
    % Plot position and force vectors, etc.
    %scatter(x,y,'*','erase','none');
    %hold on;
    %quiver(x,y,fx,fy);
    
    % Plot reference circle
    refX = all(1:Npts,2);
    refY = all(1:Npts,3);
    %scatter(refX,refY,'k','erase','none');
    %axis([-1.1 1.1 -1.1 1.1]); axis square; grid;
    %quiver(x,y,fx,fy);
    %hold off;
    
    % Find max radius and min radius and report
    Rmax(i) = max(sqrt(x.^2+y.^2));
    %Rmin = min(sqrt(x.^2+y.^2));
    
    %Frames(:,i) = getframe;
    %pause(0.25);
end

plot(1:Nframes,Rmax);
%movie(Frames,1);
