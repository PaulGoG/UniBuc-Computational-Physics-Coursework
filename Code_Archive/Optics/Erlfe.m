% Programare Ocular Erfle
% Distantele sunt in mm
clear;
clc;

function retval = translatie (d,n)
    retval = [1, 0 ; d/n , 1];
  endfunction
  
function retval = refractie (nob, nimag, R)
    retval = [1, (nob-nimag)/R ; 0 , 1];
  endfunction
  
  % Introducem datele in ordinea de pe desen, stanga - > dreapta
  % 1 dioptru desparte 2 medii => nr. dioptrii + 1 = nr. medii
  % 1 distanta desparte 2 dioptrii => nr. distante + 1 = nr. dioptrii
  
  R = [-62.89, 42.74, -42.74, 76.92, -81.30, 36.36, -36.36, 32.26, Inf] ;
  n = [1, 1.720 ,1.638 ,1 ,1.638 ,1 ,1.638 ,1.649 ,1.638, 1];
  d = [5.4 ,12.0 ,0.5, 8.0 ,0.5 ,11.5 ,1.8 ,5.5];

  Nr_refractii = length(R);
  S = [1, 0; 0, 1];
  z1 = 10;
  
    for i = 1:Nr_refractii
      
     % inmultim matricile de la coada la cap 
     indice = Nr_refractii - i + 1;
     S = S * refractie(n(indice), n(indice+1), R(indice));
     if(indice>1)
     S = S * translatie(d(indice-1), n(indice));
     endif;
     
    end;
    
    %Plane focale
    zf1 = - S(1,1)/S(1,2)
    zf2 = - S(2,2)/S(1,2)
    
    %Plane principale
    zH1 = (1- S(1,1))/S(1,2)
    zH2 = (1- S(2,2))/S(1,2)
    
    %Pozitia imaginii -> plane conjugate
    z2 = - (z1*S(2,2)+S(2,1))/(S(1,1)+z1*S(1,2))
    
    %Convergenta sistemului
    f = -1/S(1,2)
    
    %Calcul interstitiu intre planele principale
    delta = sum(d) - abs(zH1) - abs(zH2)
    
    %Marirea transversala y2 = (z2*S(1,2) + S(2,2)) * y1
      
     
     