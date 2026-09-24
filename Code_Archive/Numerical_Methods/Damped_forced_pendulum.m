% Euler-Cromer cu factor de dampare & oscilatie fortata

clear;
clc;

length = 9.8; %lungime pendul in metri
g = 9.8; 
q = 0.5; %factor dampare
F_Drive = 1.2; %
Omega_D = 2/3;
npoints = 1500;
dt = 0.04;
omega = zeros(npoints,1);
theta = zeros(npoints,1);
time = zeros(npoints,1);
theta(1) = 0.2; %cond. initiala deplasare
for step = 1:npoints-1
  omega(step+1) = omega(step) + (-(g/length)*sin(theta(step)) - q*omega(step)+F_Drive*sin(Omega_D*time(step)))*dt;
 
  tempo_theta = theta(step)+omega(step+1)*dt;
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