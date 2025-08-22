% This function creates a mixed Gaussian extensibility profile
function y = ext(beta, x)
    C1=abs(beta(1));
    C2=abs(beta(2));
    %C3=abs(beta(3));
    a1=beta(3);
    a2=beta(4);
    %a3=beta(6);
    %y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2)+C3*exp(-(x/a3).^2);
    y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2);
end

% Parameters for ext
%beta0 = [0.002,4,10,0.0075,0.001,0.0005]';
%beta0 = [0.002,4,0.0075,0.001]';
beta0 = [0.00375,0.65,0.005,0.002]';

% Visualize ext and Young's Modulus
xlim=0.01;

x=-xlim:0.0002:xlim;
y = ext(beta0, x);
z = 5000 ./ ext(beta0, x);  % Young's Modulus as function of ext
%ym = ((beta0(1) + beta0(2) - ext(beta0, x))* 1e6) + 100000;
figure, plot(x, y)
figure,plot(x, z)
%figure,plot(x, ym)

% OPTIONAL Code for a Gaussian mix in two variables
%{
% Parameters for Gaussian mix for Gaussian of two variables
A1=0.15;
A2=7;
w1=0.07;
w2=0.01;
xlim=1;
ylim=1;

% Create Gaussian mix
x=[-xlim:0.02:xlim];
y=[-ylim:0.02:ylim];
g1 = A1*exp(-(x/w1).^2);
g2 = A2*exp(-(x/w2).^2);
Gmix = g1 + g2;

% Visualize Gmix
figure,
plot(x,Gmix)

% Create Gaussian of two variables
[Xg,Yg]=meshgrid(x,y);

g1x = A1*exp(-(Xg/w1).^2);
g2x = A2*exp(-(Xg/w2).^2);
Gmixx = g1x + g2x;

g1y = A1*exp(-(Yg/w1).^2);
g2y = A2*exp(-(Yg/w2).^2);
Gmixy = g1y + g2y;

G=Gmixx.*Gmixy;

% Visualize Gaussian of two variables
figure,
surf(Xg,Yg,G)
xlabel('x')
ylabel('y')
%}