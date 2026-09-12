clear;
clc;

A=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC1.txt',' ',1,0);
B=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC2.txt',' ',1,0);
C=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC3.txt',' ',1,0);
D=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC4.txt',' ',1,0);
E=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC5.txt',' ',1,0);
F=dlmread('C:\Users\Gogu\Desktop\IOFolder\TimpiUC6.txt',' ',1,0);

subplot(2,3,1);
hist(A(:,1));
subplot(2,3,2);
hist(B(:,1));
subplot(2,3,3);
hist(C(:,1));
subplot(2,3,4);
hist(D(:,1));
subplot(2,3,5);
hist(E(:,1));
subplot(2,3,6);
hist(F(:,1));

figure

subplot(2,3,1);
hist(A(:,2));
subplot(2,3,2);
hist(B(:,2));
subplot(2,3,3);
hist(C(:,2));
subplot(2,3,4);
hist(D(:,2));
subplot(2,3,5);
hist(E(:,2));
subplot(2,3,6);
hist(F(:,2));