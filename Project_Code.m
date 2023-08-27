clear
close all
clc

di=0.03; %inner diameter of disc in m
do=0.6; %outer diameter of rotor in m
ds=0.03; %diameter of shaft in mm
E=2*10^11; %Young's Modulus in N/m2
rho=7800; %density of material in kg/m3;
t=0.05; %thickness of rotor in m;
L=2; %length of shaft in m;
g=9.81; %acceleration due to gravity

mr=rho*(pi/4)*(do^2-di^2)*t;
ms=rho*(pi/4)*(di^2)*L;
Is=(pi/64)*ds^4; %Area MOI of shaft in Nm2
Ird=(mr/4)*(do/2)^2; %Diametrical MOI of rotor
Irp=Ird/2; %Polar MOI of rotor

k1=48*E*Is/(L^3); %Bending stiffness of shaft due to weight of rotor
k2=24*E*Is/L; %Bending stiffness of shaft due to Moment of rotor on shaft
w=linspace(0,1000,1000); %omega in rpm

%CASE-1
delta=(mr*g*(L^3))/(48*E*Is);
wn_c1=sqrt(g/delta)


%CASE-3
M=zeros(4,4);
M(1,1)=mr;
M(2,2)=mr;
M(3,3)=Ird;
M(4,4)=Ird;

K=zeros(4,4);
K(1,1)=k1;
K(2,2)=k1;
K(3,3)=k2;
K(4,4)=k2;

A=inv(M)*K
[ev,lambda] = eig(K,M) %ev is the  eigen vectors, lambda is wn^2

wn1=sqrt(A(1,1)) %First Natural Frequency in rad/sec
fn1=wn1/(2*pi) %First Natural Frequency in Hertz

wn2=sqrt(A(3,3)) %Second Natural Frequency in rad/sec
fn2=wn2/(2*pi) %Second Natural Frequency in Hertz

c=Irp*(2*pi*w/60); %Varying Damping constant due to varying omega
zeta=c/(2*Ird*wn2);

wd=wn2*(sqrt(1-zeta.^2));

plot(w,wd)
xlabel('w','FontSize',15,'FontWeight','bold')
ylabel('wd','FontSize',15,'FontWeight','bold')
title('Frequency of Spin vs Damped Natural Frequency','FontSize',25)


%CASE-2
eps=0.075; %eccentricity
r=eps./((wn1./(2*pi*w/60)).^2-1);
fig=figure(2);
plot(w,r)
xlabel('w','FontSize',15,'FontWeight','bold')
ylabel('r','FontSize',15,'FontWeight','bold')
title('Spin Velocity vs Eccentricity','FontSize',25)

