% Pendulul matematic - metoda Euler

clear;
clc;

length = 1; %lungime pendul in metri
g = 9.8; 
npoints = 250;
dt = 0.04;
omega = zeros(npoints,1);
theta = zeros(npoints,1);
time = zeros(npoints,1);
theta(1) = 0.2; %cond. initiala deplasare
for step = 1:npoints-1
  omega(step+1) = omega(step) - (g/length)*theta(step)*dt;
  theta(step+1) = theta(step) + omega(step)*dt;
%   theta(step+1) = theta(step) + omega(step+1)*dt; % Euler_Cromer
  time(step+1) = time(step) + dt;
  end

  subplot(1,2,1);
  plot(time,theta,'r');
  xlabel('timp(s)');
  ylabel('theta(rad)');
  
  subplot(1,2,2);
  plot(time,omega,'b');
  xlabel('timp(s)');
  ylabel('omega');