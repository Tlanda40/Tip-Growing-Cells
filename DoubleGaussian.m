function y=DoubleGaussian(beta,x);
%Double Gaussian

x1=beta(1);
x2=beta(2);
C1=beta(3);
C2=beta(4);
a1=beta(5);
a2=beta(6);

y=C1*exp(-((x-x1)/a1).^2)+C2*exp(-((x-x2)/a2).^2);