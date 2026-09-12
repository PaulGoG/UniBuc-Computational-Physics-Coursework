clear;
clc;
%D=dlmread('TimpiUCRossler.txt',' ',1,0);
C=dlmread('DistantaLorentz.txt',' ',1,0);
D=dlmread('TimpiUCLorentz.txt',' ',1,0);
E=dlmread('DateExtreme.txt',' ',1,0);
%semilogy(C(:,2),C(:,1));
%figure
%hold on
%scatter(E(:,2),E(:,1));
%figure
edges=[0:0.05:0.8];
hist(D(:,1),edges);
figure
hist(D(:,2),edges);