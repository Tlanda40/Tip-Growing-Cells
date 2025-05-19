function y=cossquared(x,S0,a);

l=length(x);

for i=1:l
    if x(i)<=pi/2*a
        y(i)=S0*cos(x(i)/a)^2;
    elseif x(i)>pi/2*a
        y(i)=0;
    end
end