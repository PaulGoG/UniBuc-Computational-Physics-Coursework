% Programare Obiectiv Tessar
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
  
  R = [16.28, -275.70, -34.57, 15.82, Inf, 19.20, -24.00];
  n = [1, 1.6116, 1, 1.6053, 1 ,1.5123, 1.6116, 1];
  d = [3.57, 1.89, 0.81, 3.25, 2.17, 3.96];

  Nr_refractii = length(R);
  S = [1, 0; 0, 1];
 
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
    
    %Calcul interstitiu intre planele principale
    delta = sum(d) - abs(zH1) - abs(zH2)
    
    %Convergenta sistemului
    f = -1/S(1,2)    
    
    
    z1 = 3;
    y1 = 1;
    %Pozitia imaginii -> plane conjugate
    z2 = - (z1*S(2,2)+S(2,1))/(S(1,1)+z1*S(1,2))
    %Marirea transversala
    y2 = (z2*S(1,2) + S(2,2)) * y1
     
    %Marimea graficului trebuie sa cuprinda toate elementele cardinale - da eroare la plotare pentru z1 = zf1
    scalex = max([abs(z1), abs(z2), abs(zf1), abs(zf2), abs(zH1), abs(zH2), sum(d)]);
    scaley = max([abs(y1), abs(y2)]);
    
    plot([-scalex*1.25 (scalex+sum(d))*1.25], [0 0], "color", "k"); % linia
    xlim([-scalex*1.25 (scalex+sum(d))*1.25]);
    ylim([-scaley*1.5 scaley*1.5]);
    hold on;
    grid on;
    
    plot(-zH1,0,'o','MarkerFaceColor','green');text(-zH1,scaley*1.4, 'H1');
    plot([-zH1 -zH1], [-scaley*2 scaley*2], "linestyle", "-.", "color", "green");
    plot(sum(d)+zH2,0,'o','MarkerFaceColor','green');text(sum(d)+zH2,scaley*1.4, 'H2');
    plot([sum(d)+zH2 sum(d)+zH2], [-scaley*2 scaley*2], "linestyle", "-.", "color", "green");
    
    plot(0,0,'o','MarkerFaceColor','red');text(0,0.2, 'V1');
    plot(sum(d),0,'o','MarkerFaceColor','red');text(sum(d),0.2, 'V2');
    plot(-zf1,0,'o','MarkerFaceColor','blue');text(-zf1,0.2, 'F1');
    plot(sum(d)+zf2,0,'o','MarkerFaceColor','blue');text(sum(d)+zf2,0.2, 'F2');
    
    plot(-z1,0,'o','MarkerFaceColor','cyan');text(-z1,-0.1, 'ob');
    plot(-z1,y1,'^','MarkerFaceColor','cyan');
    plot([-z1 -z1], [0 y1], "linestyle", "-", "color", "cyan");
    plot(z2+sum(d),0,'o','MarkerFaceColor','magenta');text(z2+sum(d),-0.1, 'img');
    plot(z2+sum(d),y2,'^','MarkerFaceColor','magenta');
    plot([z2+sum(d) z2+sum(d)], [0 y2], "linestyle", "--", "color", "magenta");
    
     
     
     