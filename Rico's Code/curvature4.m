function [K,phi] = curvature4(X,Y)
%
%Segment lengths must be equal.
%phi is the angle from the x axis
%pphi are the internodal angles
%need to make sure normals point outwards
%This particular program is written specifically for model4.  curvature2
%and curvature3 are more general versions

[sX,temp]=size(X);

dX=diff(X);
dY=diff(Y);
Lxy=sqrt(dX.^2+dY.^2);

NNx=zeros(sX-1,1);
NNy=zeros(sX-1,1);

Xt=zeros(sX-1,1);
Yt=zeros(sX-1,1);
Lt=zeros(sX-1,1);

%Simple calculation of normals.  This only works is segment lengths are
%equal?

NNx=dX./Lxy;
NNy=dY./Lxy;

%Calculate Curvatures

pphi=zeros(sX-1,1);

for n=1:sX-1
    if NNx(n)<0 & NNy(n)<0
        pphi(n)=pi/2+atan(-NNx(n)/-NNy(n));
    elseif NNx(n)>0 & NNy(n)<0
        pphi(n)=atan(-NNy(n)/NNx(n));
    elseif NNx(n)>0 & NNy(n)>0
        pphi(n)=-atan(NNy(n)/NNx(n));
    elseif NNx(n)<0 & NNy(n)>0
        pphi(n)=-pi/2-atan(-NNx(n)/NNy(n));
    end
end

dphi=diff(pphi);

K=zeros(sX,1);

K(2:sX-1,1)=dphi./((Lxy(1:sX-2)+Lxy(2:sX-1))./2);
%K(1)=(pphi(1)-pi/2)./(Lxy(1)/2);
K(1)=K(2);
K(sX)=K(sX-1); %May want to change this

phi=zeros(sX,1);

phi(2:sX-1)=pphi(1:sX-2)+dphi/2;
phi(1)=phi(2);
phi(sX)=phi(sX-1);
