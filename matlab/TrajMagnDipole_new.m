close all;
clear all;

global B0; global q; global m; global Re; global dt; 

% parameters
e = 1.602176565e-19; % Elementary charge (Coulomb)
m_pr = 1.672621777e-27; % Proton mass (kg)
m_el = 9.10938291e-31; % Electron mass (kg)
c = 299792458; % speed of light (m/s)
% Earth parameters
B0 = 3.07e-5 % Tesla
Re = 6378137 % meter (Earth radius)

% use electron/proton
m = m_pr % we choose proton, change to m_el when simulating electron
q = e % positive charge, change to -e when simulating electron

% Trajectory of a proton with 10MeV kinetic energy in dipole field
K = 1e7    % kinetic energy in eV
K = K*e;   % convert to Joule
% Find corresponding speed:
v_mod = c/sqrt(1+(m*c^2)/K) 

% initial position: equatorial plane 4Re from Earth
%x0 = 4*Re*rand(); 
%y0 = 0.0; 
%z0 = 0.0;
R1=1.25*Re;% + 5.0*Re;
R2=3*Re;% + 5.0 * Re;
rho=rand(3);
nu=(1-2*rho(2));
x0 = (R1^3+(R2^3-R1^3)*rho(1))^(1/3)*(1-nu^2)^(1/2)*cos(2*pi*rho(3));
y0 = (R1^3+(R2^3-R1^3)*rho(1))^(1/3)*(1-nu^2)^(1/2)*sin(2*pi*rho(3));
z0 = (R1^3+(R2^3-R1^3)*rho(1))^(1/3)*nu;
pitch_angle = 30.0; % * rand() %  initial angle between velocity and mag.field (degrees)

% initial velocity
u0 = v_mod;
v0 = v_mod; %*sin(pitch_angle*pi/180); %*rand();
w0 = v_mod; %*cos(pitch_angle*pi/180); %*rand();

dt = 0.01;
tfin = 50; % in s
time = 0:dt:tfin;
%[x_sol,t] = lsode('Newton_Lorenz',[x0 ;y0; z0; u0; v0; w0],time); % solve equation of motion *** Octave ***
%[t,x_sol] = ode45('Newton_Lorenz',time,[x0 ;y0; z0; u0; v0; w0]) % use this for *** Matlab ****

x_sol = zeros(tfin/dt+1,6);
x_sol(1,1) = x0;
x_sol(1,2) = y0;
x_sol(1,3) = z0;
x_sol(1,4) = u0;
x_sol(1,5) = v0;
x_sol(1,6) = w0;

qom = q/m;

for i = 1:1:tfin/dt
	% calculate 3 Cartesian components of the magnetic field
	fac1 =  -B0*Re^3/(x_sol(i,1)^2 + x_sol(i,2)^2 + x_sol(i,3)^2)^2.5;
	Bx = 3*x_sol(i,1)*x_sol(i,3)*fac1;
	By = 3*x_sol(i,2)*x_sol(i,3)*fac1;
	Bz = (2*x_sol(i,3)^2 -x_sol(i,1)^2- x_sol(i,2)^2)*fac1;
	Ex = 0.0;
	Ey = 0.0;
	Ez = 0.0;

%    B_sq=Bx*Bx+By*By+Bz*Bz;
%    f1=tan(qom*(dt/2)*sqrt(B_sq))/sqrt(B_sq);
%    f2=2*f1/(1+f1.^2*B_sq);
%
%    vxstar2 = x_sol(i,4)+qom*Ex*0.5*dt;
%    vystar2 = x_sol(i,5)+qom*Ey*0.5*dt;
%    vzstar2 = x_sol(i,6)+qom*Ez*0.5*dt;
%
%    vxpr = vxstar2+f1*(vystar2*Bz-By*vzstar2);
%	vypr = vystar2+f1*(-vxstar2*Bz+vzstar2*Bx);
%    vzpr = vzstar2+f1*(vxstar2*By-vystar2*Bx);
%	vxplus = vxstar2+f2*(vypr*Bz-vzpr*By);
%	vyplus = vystar2+f2*(-vxpr*Bz+vzpr*Bx);
%	vzplus = vzstar2+f2*(vxpr*By-vypr*Bx);       
%	
%    x_sol(i+1,4) = vxplus+qom*0.5*dt*Ex;
%    x_sol(i+1,5) = vyplus+qom*0.5*dt*Ey;
%	x_sol(i+1,6) = vzplus+qom*0.5*dt*Ez;

	k = qom * 0.5 * dt;
	tx = k * Bx;
	ty = k * By;
	tz = k * Bz;
    t_mag2 = tx * tx + ty * ty + tz * tz;

	sx = (2.0 * tx) / (1.0 + t_mag2);
	sy = (2.0 * ty) / (1.0 + t_mag2);
	sz = (2.0 * tz) / (1.0 + t_mag2);

	v_minusx = x_sol(i,4) + k * Ex;
	v_minusy = x_sol(i,5) + k * Ey;
	v_minusz = x_sol(i,6) + k * Ez;

	v_primex = v_minusx + (v_minusy * tz - v_minusz * ty);
	v_primey = v_minusy + (v_minusz * tx - v_minusx * tz);
	v_primez = v_minusz + (v_minusx * ty - v_minusy * tx);

	v_plusx = v_minusx + (v_primey * sz - v_primez * sy);
	v_plusy = v_minusy + (v_primez * sx - v_primex * sz);
	v_plusz = v_minusz + (v_primex * sy - v_primey * sx);

	x_sol(i+1,4) = v_plusx + k * Ex;
	x_sol(i+1,5) = v_plusy + k * Ey;
	x_sol(i+1,6) = v_plusz + k * Ez;

	x_sol(i+1,1) = x_sol(i,1) + x_sol(i+1,4) * dt; % dx/dt = u
	x_sol(i+1,2) = x_sol(i,2) + x_sol(i+1,5) * dt; % dy/dt = v
	x_sol(i+1,3) = x_sol(i,3) + x_sol(i+1,6) * dt; % dz/dt = w

end

% plots
plot3(x_sol(:,1)/Re,x_sol(:,2)/Re,x_sol(:,3)/Re,'r'); % plot trajectory in 3D
xlabel('x[Re]')
ylabel('y[Re]')
zlabel('z[Re]')
title('Proton 10 MeV starting from x = 4Re, y = 0, z = 0')
hold on;
sphere(40);
hold off;

%plot(time,0.5*m*(x_sol(:,4).^2 + x_sol(:,5).^2 + x_sol(:,6).^2))

