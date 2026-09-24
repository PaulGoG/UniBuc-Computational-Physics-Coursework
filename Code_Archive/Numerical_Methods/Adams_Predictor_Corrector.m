% Adams predictor-corector in 4 pasi

clear;
clc;

% y' = f(x,y) cu y(x0) = y0
function adams4(dyfun, xspan, y0, h)
  x = xspan(1):h:xspan(2);
  % Din calcul cu functia RK obtinem 4 valori initiale
  % Folosim functia runge definita mai jos
  [xx,yy] = runge(dyfun,[x(1),x(4)],y0,h);
  y(1) = yy(1);
  y(2) = yy(2);
  y(3) = yy(3);
  y(4) = yy(4);
  
  for n=4:(length(x)-1)
    % Facem iteratii cu formula teoretica
    % feval evalueaza functia primita prin macro-ul @f in adams4 in punctele date
    p = y(n) + h/24*(55*feval(dyfun,x(n),y(n))-59*feval(dyfun,x(n-1),y(n-1))+37*feval(dyfun,x(n-2),y(n-2))-9*feval(dyfun,x(n-3),y(n-3)));
    y(n+1) = y(n) + h/24*(9*feval(dyfun,x(n+1),p)+19*feval(dyfun,x(n),y(n))-5*feval(dyfun,x(n-1),y(n-1))+feval(dyfun,x(n-2),y(n-2)));
    %display(n);
    %display(x(n));
    %display(y(n));
  endfor;
  % Comanda ' transpune vectorii din coloana-linie & vice-versa
  x = x';
  y = y';
  figure(1)
  plot(x,y,'-');
  grid on;
  xlabel('x');
  ylabel('y');
endfunction;

function res = f(x,y)
  % Functia din dreapta ODE
  res = -y+2*sin(x);
endfunction;

function [x y] = runge(dyfun,xspan,y0,h)
  % Functia va returna un vetor de vectori de lungime variabila in functie de xspan si marimea pasului
  
  x = xspan(1):h:xspan(2);
  y(1) = y0;
  for n=1:(length(x)-1)
    k1 = feval(dyfun,x(n),y(n));
    k2 = feval(dyfun,x(n)+h/2,y(n)+h*k1/2);
    k3 = feval(dyfun,x(n)+h/2,y(n)+h*k2/2);
    k4 = feval(dyfun,x(n+1),y(n)+h*k3);
    y(n+1) = y(n)+h*(k1+2*k2+2*k3+k4)/6;
  endfor;
  x = x';
  y = y';
endfunction;

adams4(@f,[0,6],0,0.1)
