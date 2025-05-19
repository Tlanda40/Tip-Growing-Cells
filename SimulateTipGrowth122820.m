%Simulates tip growth for a defined extensibility profile
%Rico Rojas
%Updated 12/28/20

%This is the simplest version of the code that simulates tip growth based
%on transverse isotropical mechancial expansion, which is equivalent to
%the model given by Equation 6 of Dumais 2004, where ass(s)=att(s), and
%ass(s)/ats(s)=nu=constant. Isotropical expansion occurs when nu=0.5.  In
%this code, ass(s) (the "extensibility profile") is called L.  The
%principal stress profiles are calculated from the apical geometry using
%Equations 1 and 2 of Dumais 2004.  From the principal stresses, the
%principal strain rates are calculated using Equation 6.  To calculate the
%displacemeent from the principal strain rates, I
%use Equations 6 and 7 from Dumais 2006 and trigonometry.  The system is
%non-dimensional, and assumes that the cell wall thickness is uniform, (and
%that the turgor pressure is constant in time).

%The simulation is initiated from an initial hemispherical geometry.

%The math is implemented directly.  Due to the discretization, the integration 
%invevitably yields a slight error in the geometry near the pole.  To
%avoid this, at each time step I slightly smooth the meridional curvature
%profile, which is then propagated to the calculation of the stresses,
%strain rates, and so on.  I've done some controls and if the smoothing is
%small enough you can avoid the instability and still get basically no
%error in the cellular geometry.  In general, the smoothing should be as
%little as possible such that you don't get an instability.  You can
%probably get away with less smoothing if you use smaller time increments.

clear
close all

profile off
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SIMULATION PARAMETERS
R=3;%Initial radius of tip
Res1=101;%Initial number of grid points to resolve apex to.
ResA=2500;%Number of grid points to pre-allocate
dt1=.05;% Initial time increment
cs=0.99;%Smoothing factor.  Prevents instabilites near the pole, caused by discretization, from exploding. (Default=0.99).
dtar=1e-2;%Simulation advances a fixed target distance, dtar, at each step by modulating time increment
nu=0.75;%Flow coupling (a.k.a. Poisson ratio)
%T=1000;%Number of time steps
T=1000;

%Define extensibility profile:
a=1;%Width of extensibility profile
L0=0.01;%Polar extensibility
%Extensibility=@cossquared;%Extensibility profile
%beta0= [0.27,10,1,0.2]';
beta0= [0.27,10,0.5,0.5]';
Extensibility=@DoubleGaussianZero;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pre-allocate matrices for variables to be tracked.  Variables are defined
%on a single half-meridian, where the first row is the pole.
r=zeros(ResA,T);%Radius
z=zeros(ResA,T);%Axial distance from pole
s=zeros(ResA,T);%Arclength from pole
phi=zeros(ResA,T);%Normal angle to surface
Ks=zeros(ResA,T);%Meridional curvature
Kt=zeros(ResA,T);%Circumferential curvature
Nx=zeros(ResA,T);%Normal vector to surface
Ny=zeros(ResA,T);%Normal vector to surface
Ss=zeros(ResA,T);%Meridional stress
St=zeros(ResA,T);%Circumferential stress
Eds=zeros(ResA,T);%Meridional strain
Edt=zeros(ResA,T);%Circumferential strain
EE=zeros(ResA,T);%Sum of strains
Vs=zeros(ResA,T);%Meridional velocity
L=zeros(ResA,T);%Extensibility
vn=zeros(ResA,T);%Normal velocity
vt=zeros(ResA,T);%Tangential velocity
ds=zeros(T,1);%Internodal distances
dt=zeros(T,1);%Time steps
Res=zeros(T,1);%Number of grid points at each timestep
Ksmid=zeros(ResA,T);%Meridional curvature at internodal points
dphi=zeros(ResA,T);%Difference in normal angle between nodes

rt=zeros(ResA,1);%Used to store displaced nodes.  That is rt(1:Res(t),1) are the displaced r(1:Res(t),t), to be resplined.
zt=zeros(ResA,1);%Same for z
st=zeros(ResA,1);%Same for s

Vsint=zeros(ResA,T);%Integrand to find meridional speed

fi=zeros(ResA,1);%Matrices used in the integration of the strain rates to find displacements of grid points
fs=zeros(ResA-1,1);
fint=zeros(ResA,1);
ff=zeros(ResA,1);

dst=zeros(T,1);%Used to store internodal distance of interpolated grid
dis=zeros(T,1);%Distance that pole advances each time step
vel=zeros(T,1);%Elongation rate
radius=zeros(T,1);%Radius of cell

%Certain variables defined on entire meridian
Ksarc=zeros(2*ResA-1,T);%Meridional curvature
Ktarc=zeros(2*ResA-1,T);%Circumferential curvature
phiarc=zeros(2*ResA-1,T);%Normal angle to surface
zsarc=zeros(2*ResA-1,T);%Axial distance from pole
rsarc=zeros(2*ResA-1,T);%Radius
sarc=zeros(2*ResA-1,T);%Arclength from pole
rtarc=zeros(ResA,1);%Used to store displaced nodes.
ztarc=zeros(ResA,1);%Used to store displaced nodes.
starc=zeros(ResA,1);%Used to store displaced nodes.
zsarcp=zeros(2*ResA-1,T);%Axial distance from pole, used to interpolate polar point
rsarcp=zeros(2*ResA-1,T);%Radius, used to interpolate polar point
ztarcp=zeros(ResA,1);%Used for re-splining polar point

sdif=zeros(2*ResA-1,T);%Arcdistance between nodes
Kerr=zeros(2*ResA-1,T);%Error obtained upon smoothing spline

%Define initial conditions
phi1=[0:pi/2/(Res1-1):pi/2]';%Normal angle to surface
r1=R*sin(phi1(1:Res1,1));%Radius
z1=R*cos(phi1(1:Res1,1));%Axial distance from pole
s1=phi1(1:Res1,1)*R;%Arclength from pole
Ks1=(1/R)*ones(1,Res1)';%Meridional curvature
Kt1=(1/R)*ones(1,Res1)';%Circumferential curvature
D1=ones(1,Res1)';%Cell wall thickness
%L1=Extensibility(s1,L0,a);%Extensibility cos functiom
L1 = Extensibility(beta0,s1);
ds1=pi/2/(Res1-1)*R;%Arclength between grid points

Res(1)=Res1;
phi(1:Res(1),1)=phi1;
r(1:Res(1),1)=r1;
z(1:Res(1),1)=z1;
s(1:Res(1),1)=s1;
Ks(1:Res(1),1)=Ks1;
Kt(1:Res(1),1)=Kt1;
D(1:Res(1),1)=D1;
L(1:Res(1),1)=L1;
ds(1)=ds1;
dt(1)=dt1;

Ksarc(1:2*Res(1)-1,1)=[flipud(Ks(1:Res(1),1))' Ks(2:Res(1),1)']';
phiarc(1:2*Res(1)-1,1)=[-flipud(phi(1:Res(1),1))' phi(2:Res(1),1)']';
zsarc(1:2*Res(1)-1,1)=[flipud(z(1:Res(1),1))' z(2:Res(1),1)']';
rsarc(1:2*Res(1)-1,1)=[-flipud(r(1:Res(1),1))' r(2:Res(1),1)']';
sarc(1:2*Res(1)-1,1)=[-flipud(s(1:Res(1),1))' s(2:Res(1),1)']';

for t=1:T-1
    
    t
    
    %In case pre-allocated matrices aren't big enough
    Eds(Res(t),t)=0;
    Edt(Res(t),t)=0;
    EE(Res(t),t)=0;
    Ss(Res(t),t)=0;
    St(Res(t),t)=0;
    
    %Calculate principal stresses using force balance
    Ss(1:Res(t),t)=1./Kt(1:Res(t),t);
    St(1:Res(t),t)=1./Kt(1:Res(t),t).*(2-Ks(1:Res(t),t)./Kt(1:Res(t),t));
    
    %Compute strain rates using the constitutive relationship
    v=ones(Res(t),1)*nu;%Define a profile for flow coupling
    Eds(1:Res(t),t)=L0*L(1:Res(t),t).*(Ss(1:Res(t),t)-v.*St(1:Res(t),t));
    Edt(1:Res(t),t)=L0*L(1:Res(t),t).*(St(1:Res(t),t)-v.*Ss(1:Res(t),t));
    EE(1:Res(t),t)=Eds(1:Res(t),t)+Edt(1:Res(t),t);
    
    %Calculate meridional speed
    Vsint(1:Res(t)-1,t)=(Eds(1:Res(t)-1,t)+Eds(2:Res(t),t)).*diff(s(1:Res(t),t))/2;
    Vs(2:Res(t),t)=cumsum(Vsint(1:Res(t)-1,t));
    
    %Calculate displacements using kinematic equations
    clear fi fs fint ff fi2 intvec intvec2
    
    fi(Res(t),1)=0;
    fs(Res(t),1)=0;
    fint(Res(t),1)=0;
    ff(Res(t),1)=0;
    ff(2:Res(t),1)=Eds(2:Res(t),t)-((Ks(2:Res(t),t)./Kt(2:Res(t),t)).*Edt(2:Res(t),t));
    fi(2:Res(t),1)=1./sin(phi(2:Res(t),t)).*ff(2:Res(t),1);
    fi2(2:Res(t)-1,1)=(fi(2:Res(t)-1,1)+fi(3:Res(t),1))*ds(t)/2;
    fint(Res(t),1)=0;
    
    intvec=flipud(fi2(2:Res(t)-1,1));
    intvec2=flipud(cumsum(intvec));
    fint(2:Res(t)-1,1)=-intvec2;
    
    vt(Res(t),t)=0;
    vn(Res(t),t)=0;
    vt(2:Res(t),t)=sin(phi(2:Res(t),t)).*fint(2:Res(t),1);
    vn(2:Res(t),t)=Edt(2:Res(t),t)./Kt(2:Res(t),t)-cos(phi(2:Res(t),t)).*fint(2:Res(t),1);%Eq16
    
    clear rt zt rtarc ztarc st starc ztarcp smid dr dz
    
    %Calcluate displacements in Euclidian coordinates
    dr=vn(2:Res(t),t)*dt(t).*sin(phi(2:Res(t),t))+vt(2:Res(t),t)*dt(t).*cos(phi(2:Res(t),t));
    dz=vn(2:Res(t),t)*dt(t).*cos(phi(2:Res(t),t))-vt(2:Res(t),t)*dt(t).*sin(phi(2:Res(t),t));
    
    %Find displaced meridians
    rt(Res(t),1)=0;
    zt(Res(t),1)=0;
    rt(2:Res(t),1)=r(2:Res(t),t)+dr;
    zt(2:Res(t),1)=z(2:Res(t),t)+dz;
    
    rtarc(2*Res(t)-1,1)=0;
    ztarc(2*Res(t)-1,1)=0;
    ztarcp(2*Res(t)-1,1)=0;
    
    rtarc(1:Res(t)-1,1)=-flipud(rt(2:Res(t),1));
    rtarc(Res(t)+1:2*Res(t)-1,1)=rt(2:Res(t),1);
    
    ztarc(1:Res(t)-1,1)=flipud(zt(2:Res(t),1));
    ztarc(Res(t)+1:2*Res(t)-1,1)=zt(2:Res(t),1);
    
    %Interpolate polar point and re-spline
    ztarcp(1:2*Res(t)-1,1)=spline([rtarc(1:Res(t)-1,1);rtarc(Res(t)+1:2*Res(t)-1,1)]',...
        [ztarc(1:Res(t)-1,1);ztarc(Res(t)+1:2*Res(t)-1,1)]',rtarc(1:2*Res(t)-1,1)');
    
    ztarc(Res(t),1)=ztarcp(Res(t),1);
    
    %Find arclengths of displaced nodes
    st(Res(t),1)=0;
    st(2:Res(t),1)=cumsum(sqrt(diff(rtarc(Res(t):2*Res(t)-1,1)).^2+diff(ztarc(Res(t):2*Res(t)-1,1)).^2));
    
    dst(t+1)=st(Res(t),1)/(Res(t)-1);
    
    starc=zeros(1,2*Res(t)-1);
    starc(1:Res(t)-1)=-flipud(st(2:Res(t),1))';
    starc(Res(t):2*Res(t)-1)=st(1:Res(t),1)';
    intarc1=(-st(Res(t),1):dst(t+1):st(Res(t),1));
    
    clear rparc zparc sp sparc
    
    %Find resolution of displaced meridian
    Res(t+1)=ceil(-starc(1,1)/ds(1));
    
    %Find new internodal arclength
    ds(t+1)=-starc(1,1)/(Res(t+1)-1);
    
    %Re-spline displaced meridian
    rsarcp(2*Res(t+1)-1,t+1)=0;
    zsarcp(2*Res(t+1)-1,t+1)=0;
    intarc=(starc(1,1):ds(t+1):-starc(1,1));    
    rsarcp(1:2*Res(t+1)-1,t+1)=spline(starc(1,1:2*Res(t)-1)',rtarc(1:2*Res(t)-1,1)',intarc)';
    zsarcp(1:2*Res(t+1)-1,t+1)=spline(starc(1,1:2*Res(t)-1)',ztarc(1:2*Res(t)-1,1)',intarc)';
    
    rsarc(2*Res(t+1)-1,t+1)=0;
    zsarc(2*Res(t+1)-1,t+1)=0;
    rsarc(1:2*Res(t+1)-1,t+1)=rsarcp(Res(t+1)-Res(t+1)+1:Res(t+1)+Res(t+1)-1,t+1);
    zsarc(1:2*Res(t+1)-1,t+1)=zsarcp(Res(t+1)-Res(t+1)+1:Res(t+1)+Res(t+1)-1,t+1);
    
    r(Res(t+1),t+1)=0;
    z(Res(t+1),t+1)=0;
    s(Res(t+1),t+1)=0;
    r(1:Res(t+1),t+1)=rsarc(Res(t+1):2*Res(t+1)-1,t+1);
    z(1:Res(t+1),t+1)=zsarc(Res(t+1):2*Res(t+1)-1,t+1);
 
    sdif(1:Res(t+1)-1,t+1)=sqrt(diff(r(1:Res(t+1),t+1)).^2+diff(z(1:Res(t+1),t+1)).^2);
    s(2:Res(t+1),t+1)=cumsum(sdif(1:Res(t+1)-1,t+1));
    
    sarc(2*Res(t+1)-1,t+1)=0;
    sarc(1:Res(t+1)-1,t+1)=-flipud(s(2:Res(t+1),t+1));
    sarc(Res(t+1):2*Res(t+1)-1,t+1)=s(1:Res(t+1),t+1);
    
    %Calculate polar displacement and elongation rate
    dis(t)=z(1,t+1)-z(1,t);
    vel(t)=dis(t)/dt(t);
    
    %Compute new geometrical properties
    Ks(Res(t+1),t+1)=0;
    Kt(Res(t+1),t+1)=0;
    phi(Res(t+1),t+1)=0;
    Ksarc(2*Res(t+1)-1,t+1)=0;
    Ktarc(2*Res(t+1)-1,t+1)=0;
    phiarc(2*Res(t+1)-1,t+1)=0;
    
    [Ksarc(1:2*Res(t+1)-1,t+1),phiarc(1:2*Res(t+1)-1,t+1)]=curvature4(rsarc(1:2*Res(t+1)-1,t+1),zsarc(1:2*Res(t+1)-1,t+1));
    
    %Smooth meridional curvature profile
    Kpre=Ksarc(1:2*Res(t+1)-1,t+1)';
    Ksarc(1:2*Res(t+1)-1,t+1)=csaps(sarc(1:2*Res(t+1)-1,t+1)',Kpre,cs,sarc(1:2*Res(t+1)-1,t+1)')';
    
    %Calculate error in curvature
    Kerr(1:2*Res(t+1)-1,t+1)=Ksarc(1:2*Res(t+1)-1,t+1)-Kpre';
    
    %Recalculate normal angle based on smoothed Ks.
    Ks(1:Res(t+1),t+1)=Ksarc(Res(t+1):2*Res(t+1)-1,t+1);
    
    Ksmid(1:Res(t+1)-1,t+1)=(Ks(1:Res(t+1)-1,t+1)+Ks(2:Res(t+1),t+1))/2;
    dphi(1:Res(t+1)-1,t+1)=Ksmid(1:Res(t+1)-1,t+1).*sdif(1:Res(t+1)-1,t+1);
    
    phi(2:Res(t+1),t+1)=cumsum(dphi(1:Res(t+1)-1,t+1));
    phiarc(1:2*Res(t+1)-1,t+1)=[-flipud(phi(1:Res(t+1),t+1));phi(2:Res(t+1),t+1)];
    
    Kt(2:Res(t+1),t+1)=sin(phi(2:Res(t+1),t+1))./r(2:Res(t+1),t+1);
    Kt(1,t+1)=Ks(1,t+1);
    Ktarc(1:Res(t+1)-1,t+1)=flipud(Kt(2:Res(t+1),t+1));
    Ktarc(Res(t+1)+1:2*Res(t+1)-1,t+1)=Kt(2:Res(t+1),t+1);
    Ktarc(Res(t+1),t+1)=Ksarc(Res(t+1),t+1);
    
    %Calculate extensibility profile for next time step
    L(Res(t+1),t+1)=0;
    L(1:Res(t+1),t+1)=Extensibility(beta0,s(1:Res(t+1),t+1));%DoubleGaussian
    %L(1:Res(t+1),t+1)=Extensibility(s(1:Res(t+1),t+1),L0,a);%%-cosfunction
    
    %Modulate time increment to achieve target polar displacement
    if dis(t)>1e-6
        dt(t+1)=dtar/dis(t)*dt(t);
    else 
        dt(t+1)=dt(t);
    end
end

toc

%Calculate total time
tt=sum(dt);
time=cumsum(dt(2:T));

%Plot results

%nth=50;%Plot every nth outline
nth=50;

f1=figure(1);
 hold on
for t=1:nth:T
    %figure
    plot(zsarc(1:2*Res(t)-1,t),rsarc(1:2*Res(t)-1,t),'-');

    axis equal
    axis off
    xlim([0 13])
    grid on
    fig2pretty
end
% axis equal
% axis off
% grid on
% fig2pretty

% f2=figure(2);
% clf
% hold on
% for t=1:nth:T
%     plot(sarc(1:2*Res(t)-1,t),Ksarc(1:2*Res(t)-1,t),'-');
% end
% plot(sarc(1:2*Res(T)-1,T),Ksarc(1:2*Res(T)-1,T),'-r','Linewidth',2)
% xlabel('Meriodonal Arclength')
% ylabel('Meridional Curvature')
% fig2pretty
% 
% f4=figure(4);
% set(f4,'NumberTitle','off','Name','Meriodional Strain rate')
% clf
% hold on
% for t=2:nth:T-1
%     plot(s(1:Res(t),t),Eds(1:Res(t),t),'-');
% end
% grid on
% xlabel('Meriodonal Arclength')
% ylabel('Meridional Strain Rate')
% fig2pretty
% 
% f5=figure(5);
% set(f5,'NumberTitle','off','Name','Elongation Rate')
% clf
% hold on
% plot(time(1:T-1), vel(1:T-1));
% grid on
% xlabel('Time')
% ylabel('Elongation Rate')
% fig2pretty
% 
% f7=figure(7);
% set(f7,'NumberTitle','off','Name','Error in Curvature')
% clf
% hold on
% plot(sarc(1:2*Res(T)-1,T),Kerr(1:2*Res(T)-1,T),'-r')
% grid on
% xlabel('Meriodonal Arclength')
% ylabel('Error in Meridional Curvature')
% fig2pretty