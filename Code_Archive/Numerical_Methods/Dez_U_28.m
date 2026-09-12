% Dezintegrarea Uraniului 238 - Metoda Euler
% Ecuatia: dN/dt = -N/tau
clear;
clc;

N_uranium_initial = 1000; %nr initial nuclee
npoints = 1000; %discretizarea intervalului temporal
dt = 1e7; %unitate de pas in ani
tau = 4.4e9; %timpul mediu de viata U238
N_uranium = zeros(npoints, 1 ); %declaram matrice de marime npoints x 1, umpluta cu 0-uri
time = zeros(npoints, 1 ); %declaram vectorul timp
N_uranium(1) = N_uranium_initial; %conditia initiala
time(1) = 0; %alta conditie initiala

  for step=1:npoints-1
    N_uranium(step+1) = N_uranium(step) - (N_uranium(step)/tau) * dt;
    time(step+1) = time(step) + dt;
  end
  
  subplot(1,3,1);
  plot(time, N_uranium,'r');
  xlabel('Timpul in ani de zile');
  ylabel('Numar nuclee');
  title('Numeric');

  t = 0:1e8:10e9;
  N_analytical = N_uranium_initial*exp(-t/tau);
  subplot(1,3,2);
  plot(t,N_analytical,'b');
  xlabel('Timpul in ani de zile');
  ylabel('Numar nuclee');
  title('Analitic');
  
  subplot(1,3,3);
  plot(time, N_uranium,'r');
  hold;
  scatter(t,N_analytical,'b');
  xlabel('Timpul in ani de zile');
  ylabel('Numar nuclee');
  title('Comparatie');