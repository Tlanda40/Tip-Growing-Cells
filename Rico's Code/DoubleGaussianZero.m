function y=DoubleGaussianZero(beta,x);
%Double Gaussian

C1=abs(beta(1));
C2=abs(beta(2));
a1=beta(3);
a2=beta(4);

y=C1*exp(-(x/a1).^2)+C2*exp(-(x/a2).^2);