clear;
clc;

A=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC.txt',' ',1,0);

pkg load optim;

[y,x]=hist(A(:,1),10);

F=@(x, p) p(1)*p(2)*x.^(-p(3))+(1-p(1))*p(4)*exp(-(x-p(5)).^2/p(6));
pin=[rand,rand,rand+1,rand,rand,rand];

B=[0,1;0,Inf;1,Inf;0,Inf;0,Inf;0,Inf];
options.bounds=B;

[f,p,cvg,iter]=leasqr(x,y,pin,F,.00001,100,ones(size(y)),.001*ones(size(pin)),"dfdp",options);

Functie=@(x) p(1)*p(2)*x.^(-p(3))+(1-p(1))*p(4)*exp(-(x-p(5)).^2/p(6));

hist(A(:,1),10);
hold on
ezplot(Functie, [0,1]);

p