%Runge-Kutta pentru sistem de ODE

clear;
clc;

function dy = f(t,y)
  dy = [y(2);
  -y(1)-2*exp(t)+1;
  -y(1)-exp(t)+1];
endfunction

function [t w] = rk4_systems(a,b,N,alpha)
  
  m = size(alpha);
  if m==1
    alpha = alpha';
  endif
  
  h = (b-a)/N;
  t(1) = a;
  w(:,1) = alpha;
  
  for i = 1:N
    k1 = h*f(t(i), w(:,i));
    k2 = h*f(t(i)+h/2, w(:,i)+0.5*k1);
    k3 = h*f(t(i)+h/2, w(:,i)+0.5*k2);
    k4 = h*f(t(i)+h,w(:,i)+k3);
    w(:,i+1) = w(:,i)+(k1+2*k2+2*k3+k4)/6;
    t(i+1) = a+i*h;
  endfor
  [t' w']
endfunction

a = 1;
b = 8;
no = 100;
init = [1 2 3];
[t w] = rk4_systems(a,b,no,init)

n = length(t);
i = 1;
contor = 0;

while i<=3*n
  contor = contor+1;
  y1sol(contor) = w(i);
  i = i+1;
  y2sol(contor) = w(i);
  i = i+1;
  y3sol(contor) = w(i);
  i = i+1;
endwhile


  subplot(1,3,1);
  plot(t,y1sol,'r');
  xlabel('x');
  ylabel('y1');
  
  subplot(1,3,2);
  plot(t,y2sol,'b');
  xlabel('x');
  ylabel('y2');
  
  subplot(1,3,3);
  plot(t,y3sol,'g');
  xlabel('x');
  ylabel('y3');
  
  figure;
  plot3(y1sol,y2sol,y3sol,'m');
  xlabel('y1');
  ylabel('y2');
  zlabel('y3');
  
  