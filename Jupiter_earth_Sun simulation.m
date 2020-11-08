# This Code simulates the motion of Jupiter and earth with respect to sun (at center). 

npoints=1000000;
dt = 0.0001; % time step in years.
M_s=2e30; % Mass of the Sun in kg
M_e1=6e24; % Mass of the Earth in kg
M_e=1*M_e1; % playing with the mass of earth
Mj_actual=1.9e27; % Actual mass of Jupiter
M_j=1*Mj_actual; % playing with the mass of jupiter
x_e_initial=1; % Initial position of Earth in AU
y_e_initial=0;
v_e_x_initial=0; % Initial velocity of Earth in AU/yr
v_e_y_initial=2*pi;
x_j_initial=5.2; % Initial position of Jupiter in AU
y_j_initial=0;
v_j_x_initial=0; % Initial velocity of Jupiter in AU/yr
v_j_y_initial= 2.7549; % This is 2*pi*5.2 AU/11.85 years = 2.75 AU/year
% Create arrays to store position and velocity of Earth
x_e=zeros(npoints,1);
y_e=zeros(npoints,1);
v_e_x=zeros(npoints,1);
v_e_y=zeros(npoints,1);
x_j=zeros(npoints,1);
y_j=zeros(npoints,1);
v_j_x=zeros(npoints,1);
v_j_y=zeros(npoints,1);
r_e=zeros(npoints,1);
r_j=zeros(npoints,1);
r_e_j=zeros(npoints,1);
% Initialise positions and velocities of Earth and Jupiter
x_e(1)=x_e_initial;
y_e(1)=y_e_initial;
v_e_x(1)=v_e_x_initial;
v_e_y(1)=v_e_y_initial;
x_j(1)=x_j_initial;
y_j(1)=y_j_initial;
v_j_x(1)=v_j_x_initial;
v_j_y(1)=v_j_y_initial;
for i = 1:npoints-1 % loop over the timesteps
% Calculate distances to Earth from Sun, Jupiter from Sun and Jupiter to Earth for current value of i
r_e(i)=sqrt(x_e(i)^2+y_e(i)^2);
r_j(i)=sqrt(x_j(i)^2+y_j(i)^2);
r_e_j(i)=sqrt((x_e(i)-x_j(i))^2 +(y_e(i)-y_j(i))^2);
v_e_x(i+1)=v_e_x(i)-4*pi^2*x_e(i)*dt/r_e(i)^3-4*pi^2*(M_j/M_s)*(x_e(i)-x_j(i))*dt/r_e_j(i)^3;
v_e_y(i+1)=v_e_y(i)-4*pi^2*y_e(i)*dt/r_e(i)^3-4*pi^2*(M_j/M_s)*(y_e(i)-y_j(i))*dt/r_e_j(i)^3;
% Compute x and y components for new velocity of Jupiter
v_j_x(i+1)=v_j_x(i)-4*pi^2*x_j(i)*dt/r_j(i)^3-4*pi^2*(M_e/M_s)*(x_j(i)-x_e(i))*dt/r_e_j(i)^3;
v_j_y(i+1)=v_j_y(i)-4*pi^2*y_j(i)*dt/r_j(i)^3-4*pi^2*(M_e/M_s)*(y_j(i)-y_e(i))*dt/r_e_j(i)^3;
% Use Euler  technique to calculate the new positions of Earth and Jupiter. 
x_e(i+1)=x_e(i)+v_e_x(i+1)*dt;
y_e(i+1)=y_e(i)+v_e_y(i+1)*dt;
x_j(i+1)=x_j(i)+v_j_x(i+1)*dt;
y_j(i+1)=y_j(i)+v_j_y(i+1)*dt;
end
plot(x_e,y_e, 'r', x_j,y_j, 'k');
axis([-7 7 -7 7]);
xlabel('x(AU)');
ylabel('y(AU)');
title('3 body simulation - Jupiter-Earth-sun');
