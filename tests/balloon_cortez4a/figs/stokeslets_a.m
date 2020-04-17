% files
file1 = 'p_1';
file2 = 'p_2';
file3 = 'p_3';
file4 = 'p_4';
infile_1 = 'pressure_exact.txt';
infile_2 = 'pressure_stokeslets.txt';
infile_3 = 'u_exact.txt';
infile_4 = 'v_exact.txt';

% draw exact pressure function along the line (x, 3/10)
x = 0:0.01:2;
p = ...
( -3*0.3*x.^2. + 0.3.^3. ) .* ( x < sqrt(91)/10) + ...
( (3*0.3*x.^2. - 0.3.^3.)./(x.^2. + 0.3.^2.).^3. ) .* ( x > sqrt(91)/10);

% draw exact u function along the line (x, 3/10)
r = sqrt(x.^2. + 0.3.^2);
t = atan2(0.3, x);
u = ...
( (3/8).*(r.^2.).*sin(2.*t) + (1/16).*(r.^4.).*sin(4.*t) - (1/4).*(r.^4.).*sin(2.*t) ) .* ( x < sqrt(91)/10 ) + ...
( (1/8).*(r.^-2.).*sin(2.*t) - (3/16).*(r.^-4.).*sin(4.*t) + (1/4).*(r.^-2.).*sin(4.*t) ) .* ( x >= sqrt(91)/10);

% draw exact v function along the line (x, 3/10)
v = ...
( (3/8).*(r.^2.).*cos(2.*t) - (1/16).*(r.^4.).*cos(4.*t) - (1/4).*(r.^4.).*cos(2.*t) ) .* ( x < sqrt(91)/10 ) + ...
( (1/8).*(r.^-2.).*cos(2.*t) + (3/16).*(r.^-4.).*cos(4.*t) - (1/4).*(r.^-2.).*cos(4.*t) ) .* ( x>= sqrt(91)/10 );

% draw computed pressure p along the line (x, 3/10) for N = 100
Z1 = load(infile_1);
X1 = Z1(:,1);
P1 = Z1(:,2);

plot(X1,P1,'go',x,p);

axis([ 0 2 -0.8 0.8 ]);
saveas(gcf,file1,'tiff');

% draw computed pressure using stokeslets along the line (x, 3/10) for N = 100
Z2 = load(infile_2);
X2 = Z2(:,1);
P2 = Z2(:,2);

plot(X2,P2,'ro',x,p);

axis([ 0 2 -0.8 0.8 ]);
saveas(gcf,file2,'tiff');

% draw computed u velocity
Z3 = load(infile_3);
X3 = Z3(:,1);
U3 = Z3(:,2);

plot(X3,U3,'go',x,u);

axis([ 0 2 0 0.14 ]);
saveas(gcf,file3,'tiff');

% draw computed v velocity
Z4 = load(infile_4);
X4 = Z4(:,1);
V4 = Z4(:,2);

plot(X4,V4,'ro',x,v);

axis([ 0 2 -0.04 0.12 ]);
saveas(gcf,file4,'tiff');

