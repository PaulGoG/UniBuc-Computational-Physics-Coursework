% Runge-Kutta 4 cu factor de dampare & oscilatie fortata

clear;
clc;

function res = frhs(omega,theta,time)
  
length = 9.8; %lungime pendul in metri
g = 9.8; 
q = 0.5; %factor dampare
F_Drive = 1.2; %
Omega_D = 2/3;

res = (-g/(length)*sin(theta)-q*omega+F_Drive*sin(Omega_D*time));
endfunction

npoints = 1500;
dt = 0.04;
omega = zeros(npoints,1);
theta = zeros(npoints,1);
time = zeros(npoints,1);
theta(1) = 0.2; %cond. initiala deplasare
omega(1) = 0;


for step = 1:npoints-1
  
  k1 = frhs(omega(step), theta(step), time(step));
  k2 = frhs(omega(step) + dt*k1/2, theta(step) + dt*k1/2, time(step) + dt/2);
  k3 = frhs(omega(step) + dt*k2/2, theta(step) + dt*k2/2, time(step) + dt/2);
  k4 = frhs(omega(step) + dt*k3, theta(step) + dt*k3, time(step) + dt);
  
  omega(step+1) = omega(step) + (1/6)*dt*(k1+k2*2+k3*2+k4);
  
  pk1 = omega(step);
  pk2 = omega(step) + dt*pk1/2.0;
  pk3 = omega(step) + dt*pk2/2.0;
  pk4 = omega(step) + dt*pk3;
  
  tempo_theta = theta(step)+ (1.0/6.0)*dt*(pk1+pk2*2+pk3*2+pk4);
  
  %theta intre - pi si pi pentru reprezentarile grafice
  
  if(tempo_theta<-pi)
  tempo_theta=tempo_theta+2*pi;
  elseif(tempo_theta>pi)
  tempo_theta=tempo_theta-2*pi;
  end;


  theta(step+1) = tempo_theta;
  time(step+1) = time(step) + dt;
  end;


  subplot(1,3,1);
  plot(theta,omega,'r');
  xlabel('theta(rad)');
  ylabel('omega(s)');
  
  subplot(1,3,2);
  plot(time,theta,'b');
  xlabel('timp(s)');
  ylabel('theta(rad)');
  
  subplot(1,3,3);
  plot(time,omega,'g');
  xlabel('timp(s)');
  ylabel('omega(s)');