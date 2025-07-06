A1=0.27;
A2=10;
w1=0.1;
w2=0.1;
xlim=1;
ylim=1;

x=[-xlim:0.02:xlim];
y=[-ylim:0.02:ylim];
g1 = A1*exp(-(x/w1).^2);
g2 = A2*exp(-(x/w2).^2);
Gmix = g1 + g2;

[Xg,Yg]=meshgrid(x,y);

g1x = A1*exp(-(Xg/w1).^2);
g2x = A2*exp(-(Xg/w2).^2);
Gmixx = g1x + g2x;

g1y = A1*exp(-(Yg/w1).^2);
g2y = A2*exp(-(Yg/w2).^2);
Gmixy = g1y + g2y;

G=Gmixx.*Gmixy;



figure,
surf(Xg,Yg,G)
xlabel('x')
ylabel('y')

figure,plot(x,Gmix), fig2pretty
