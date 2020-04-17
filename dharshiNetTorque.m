% Torque computation using Dharshi's energy functional for single triangle
% expectation: zero net torque

clear all;

%% constants
lambda = 1.0;
mu = 1.0;

%% initial configuration
s = rand(3,2);
s0 = s(1,:);
s1 = s(2,:);
s2 = s(3,:);

%% deformed configuration
e = 0.02;
%x = (1+e)*s; % linear deformation
x(:,1) = s(:,1) + e*s(:,1).*s(:,1); % nonlinear deformation
x(:,2) = s(:,2); % nonlinear deformation
x0 = x(1,:);
x1 = x(2,:);
x2 = x(3,:);

%% Compute elastic force
forces = zeros(3,2);

% Vref term (3)
det = (s1(:,1)-s0(:,1))*(s2(:,2)-s0(:,2)) - (s2(:,1)-s0(:,1))*(s1(:,2)-s0(:,2));
dA = 0.5*abs(det);

% invS term (2)
t2_0_x = zeros(2,2);
t2_0_x(1,1) = (1/det)*(s1(:,2)-s2(:,2));
t2_0_x(1,2) = (1/det)*(s2(:,1)-s1(:,1));
t2_0_x(2,1) = 0;
t2_0_x(2,2) = 0;

t2_0_y = zeros(2,2);
t2_0_y(1,1) = 0;
t2_0_y(1,2) = 0;
t2_0_y(2,1) = (1/det)*(s1(:,2)-s2(:,2));
t2_0_y(2,2) = (1/det)*(s2(:,1)-s1(:,1));

t2_1_x = zeros(2,2);
t2_1_x(1,1) = (1/det)*(s2(:,2)-s0(:,2));
t2_1_x(1,2) = (1/det)*(s0(:,1)-s2(:,1));
t2_1_x(2,1) = 0;
t2_1_x(2,2) = 0;

t2_1_y = zeros(2,2);
t2_1_y(1,1) = 0;
t2_1_y(1,2) = 0;
t2_1_y(2,1) = (1/det)*(s2(:,2)-s0(:,2));
t2_1_y(2,2) = (1/det)*(s0(:,1)-s2(:,1));

t2_2_x = zeros(2,2);
t2_2_x(1,1) = (1/det)*(s0(:,2)-s1(:,2));
t2_2_x(1,2) = (1/det)*(s1(:,1)-s0(:,1));
t2_2_x(2,1) = 0;
t2_2_x(2,2) = 0;

t2_2_y = zeros(2,2);
t2_2_y(1,1) = 0;
t2_2_y(1,2) = 0;
t2_2_y(2,1) = (1/det)*(s0(:,2)-s1(:,2));
t2_2_y(2,2) = (1/det)*(s1(:,1)-s0(:,1));

% deformation gradient tensor, A
a = zeros(2,2);
a(1,1) = (1/det)*( (x1(:,1)-x0(:,1))*(s2(:,2)-s0(:,2)) + (x2(:,1)-x0(:,1))* ...
    (s0(:,2)-s1(:,2)) );
a(1,2) = (1/det)*( (x1(:,1)-x0(:,1))*(s0(:,1)-s2(:,1)) + (x2(:,1)-x0(:,1))* ...
    (s1(:,1)-s0(:,1)) );
a(2,1) = (1/det)*( (x1(:,2)-x0(:,2))*(s2(:,2)-s0(:,2)) + (x2(:,2)-x0(:,2))* ...
    (s0(:,2)-s1(:,2)) );
a(2,2) = (1/det)*( (x1(:,2)-x0(:,2))*(s0(:,1)-s2(:,1)) + (x2(:,2)-x0(:,2))* ...
    (s1(:,1)-s0(:,1)) );

% first Piola-Kirchhoff term (1)
t1 = zeros(2,2);
t1(1,1) = (a(1,1)-1)*(lambda+2*mu) + lambda*(a(2,2)-1);
t1(1,2) = mu*(a(1,2)+a(2,1));
t1(2,1) = t1(1,2);
t1(2,2) = (a(2,2)-1)*(lambda+2*mu) + lambda*(a(1,1)-1);

%
forces(1,1) = -dA*(t1(1,1)*t2_0_x(1,1) + t1(1,2)*t2_0_x(1,2) + t1(2,1)*t2_0_x(2,1) + t1(2,2)*t2_0_x(2,2));
forces(1,2) = -dA*(t1(1,1)*t2_0_y(1,1) + t1(1,2)*t2_0_y(1,2) + t1(2,1)*t2_0_y(2,1) + t1(2,2)*t2_0_y(2,2));

forces(2,1) = -dA*(t1(1,1)*t2_1_x(1,1) + t1(1,2)*t2_1_x(1,2) + t1(2,1)*t2_1_x(2,1) + t1(2,2)*t2_1_x(2,2));
forces(2,2) = -dA*(t1(1,1)*t2_1_y(1,1) + t1(1,2)*t2_1_y(1,2) + t1(2,1)*t2_1_y(2,1) + t1(2,2)*t2_1_y(2,2));

forces(3,1) = -dA*(t1(1,1)*t2_2_x(1,1) + t1(1,2)*t2_2_x(1,2) + t1(2,1)*t2_2_x(2,1) + t1(2,2)*t2_2_x(2,2));
forces(3,2) = -dA*(t1(1,1)*t2_2_y(1,1) + t1(1,2)*t2_2_y(1,2) + t1(2,1)*t2_2_y(2,1) + t1(2,2)*t2_2_y(2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute torque using reference coordinates
Tref = s0(:,1)*forces(1,2)-s0(:,2)*forces(1,1) + s1(:,1)*forces(2,2)-s1(:,2)*forces(2,1) + ...
    s2(:,1)*forces(3,2)-s2(:,2)*forces(3,1);

%% Compute torque using deformed coordinates
Tdef = x0(:,1)*forces(1,2)-x0(:,2)*forces(1,1) + x1(:,1)*forces(2,2)-x1(:,2)*forces(2,1) + ...
    x2(:,1)*forces(3,2)-x2(:,2)*forces(3,1);

fprintf('torque reference: %.16f, torque deformed: %.16f\n',Tref,Tdef);