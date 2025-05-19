%DoubleGaussianPeak.m
%
clear, close all

A1=0.27;
A2=10;
w1=0.5;
w2=0.5;
xlim=20;
ylim=5;

x=[-xlim:0.15:xlim];
y=[-ylim:0.15:ylim];
g1 = A1*exp(-(x/w1).^2);
g2 = A2*exp(-(x/w2).^2);
Gmix = g1 + g2;

[Xg,Yg]=meshgrid(x,y);

G1=A1*exp(-(Xg/w1).^2);
G2=A2*exp(-(Yg/w2).^2);

G=G1.*G2;



figure,
surf(Xg,Yg,G)
xlabel('g1')
ylabel('g2')

figure,plot(x,Gmix), fig2pretty
