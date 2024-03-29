%% This code is for the 2D advection-diffusion equation or Burgers' Equ.
%% A Matlab code written by and developed by HUSSEIN A. H. Muhammed Nov.-Dec. 2022.
%% B.Sc.H AND M.Sc. (Honuors).
%% This code generates a video file that simulates the simplest form of Navier-Stocks' Equ. in Cartesian coordinate system
%% Using Advection and Diffusion forces. Thanks is to Dr. Haroon Stephen (UNLV-Las Vegas, USA).
% Lastest update 04. Mar. 2024

clear;

%%Equation
%% Ct = K*[ Cxx + Cyy ] - u [Cx + Cy] + qc
%% where C is the scalar quantity being transported,
%% K=diffusion coefficient;  u  is the velocity
%% components in the x and y directions, which is equal. 

% %% Prepare the movie file
%     vidObj = VideoWriter('BE-2d.avi');
%     open(vidObj);

%% Domain
% Space
Lx=7;
Ly=10;
dx=0.2; dy=dx;
nx=fix(Lx/dx);
ny=fix(Ly/dy); 
X=linspace(0, Lx, nx);
Y=linspace(0, Ly, ny);
[x,y]=meshgrid(X,Y);
x=x'; y=y';

%% Propagation Time
T=8;

%% Field Arrays
% Variables
C=zeros(nx,ny);

%% Parameters
u=zeros(nx,ny);
v=zeros(nx,ny);
K=zeros(nx,ny);


%% Intial Conditions
t=0;
C(:)=0;
u(:)=1.5;
v=0.3+0.01*(y-Ly)+cos(4*pi*x/Lx);

CFL=0.5;
dt=CFL*min(dx./abs(u(:)) + dy./abs(v(:)));


%% Time stepping while Loop
while (t < T)
    % Zero B.C.
    C(:,[1 end]) = C(:,[2 end-1]);
    C([1 end],:) = C([2 end-1],:);

    % Oscilation Source
    if t < 2
        C(4,12)=C(4,12)+dt+50;
    end

    % Conventional Space-Time Descritized Solution
    t=t+dt;
    Cn=C;
    for i=2:nx-1,  for j=2:ny-1 
        % Advection term
        A = u(i,j)*(Cn(i+1,j) - Cn(i-1,j))/(2*dx)  ...
          + v(i,j)*(Cn(i,j+1) - Cn(i,j-1))/(2*dy);
        % Diffusion term
        D = K(i,j)*(Cn(i+1,j)-2*Cn(i,j)+Cn(i-1,j))/dx^2 ...
          + K(i,j)*(Cn(i,j+1)-2*Cn(i,j)+Cn(i,j-1))/dy^2;
        % Euler's Method
        C(i,j) = Cn(i,j) + dt*(-A + D);            % avoid negative velocity values
        if C(i,j) < 0; C(i,j)=0; end
    end,  end

%% 2D Visualize at selected steps
% clf;
% imagesc(X, Y, C');
% colorbar; 
% c = colorbar;
% c.Label.String = 'Velocity Field (unit/sec.)';
% hold on
% quiver(x,y,u,v);
% hold off;
% set(gca, 'ydir', 'norm');
% title(sprintf('propagation time = %.2f' , t));
% axis([0 Lx 0 Ly]);
% xlabel('horizental axis');
% ylabel('vertical axis');
% zlabel('velocity field axis');
% legend;
% shg; pause(0.1);

%% 3D Visualize at selected steps
clf;
surf(X, Y, C');  % Create 3D surface plot
colorbar;
c = colorbar;
c.Label.String = 'Velocity Field (unit/sec.)';
hold on;
quiver3(x, y, zeros(size(u)), u, v, zeros(size(v)), 'k');  % Plot velocity vectors on the surface
hold off;
title(sprintf('Propagation Time = %.2f', t));
xlabel('Horizontal Axis');
ylabel('Vertical Axis');
zlabel('Velocity Field Axis');
axis([0 Lx 0 Ly min(min(C)) max(max(C))]);  % Set axis limits
legend('Velocity Vectors');
shg;
pause(0.1);


 % % Write each frame to the file
 %       currFrame = getframe(gcf);
 %       writeVideo(vidObj,currFrame);
end


% %% Close the file
% close(vidObj);
