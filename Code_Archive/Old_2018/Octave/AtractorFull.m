clear;
clc;

A=dlmread('C:\Users\Gogu\Desktop\IOFolder\Atractor1.txt',' ',1,0);
B=dlmread('C:\Users\Gogu\Desktop\IOFolder\Atractor2.txt',' ',1,0);
C=dlmread('C:\Users\Gogu\Desktop\IOFolder\Distanta.txt',' ',1,0);
D=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC.txt',' ',1,0);

plot3(A(:,1),A(:,2),A(:,3));
figure
plot3(B(:,1),B(:,2),B(:,3));
figure
plot(C(:,2),C(:,1));
figure
hist(D(:,1));
figure
hist(D(:,2));