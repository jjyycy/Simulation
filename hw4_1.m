clc
clear
% a
n=10000;
S0=100;
K=100;
r=0.05;
T=1;
v0square=0.04;
N=50;
alpha=0.1;
phi=0.1; % or 1.0
delta=T/N;

S=S0*ones(n,1);
vsquare=v0square*ones(n,1);

for i=1:N
    z=randn(n,2);
    S=S+r*delta*S+sqrt(delta)*sqrt(vsquare).*S.*z(:,1);
    vsquare=vsquare+alpha*delta.*vsquare+phi*sqrt(delta).*vsquare.*z(:,2);
end

C_a=exp(-r*T).*max(S-K,0);
Cbar_a=mean(C_a)
Cstderr_a=std(C_a)/sqrt(n)

% b
vsquare2=zeros(N+1,n);
vsquare2(1,:)=v0square*ones(1,n);

for i=1:N
    z=randn(1,n);
    vsquare2(i+1,:)=vsquare2(i,:)+alpha*delta.*vsquare2(i,:)+phi*sqrt(delta).*vsquare2(i,:).*z;
end

vol=zeros(1,n);

for i=2:(N+1)
    vol=vol+vsquare2(i,:);
end
vol=sqrt(vol/N);

% Black-Scholes formula:
dplus=1./(vol.*sqrt(T)).*(log(S0/K)+(r+1/2.*vol.^2)*T);
dminus=1./(vol.*sqrt(T)).*(log(S0/K)+(r-1/2.*vol.^2)*T);
C_b=S0*normcdf(dplus)-normcdf(dminus)*K*exp(-r*T);
Cbar_b=mean(C_b)
Cstderr_b=std(C_b)/sqrt(n)