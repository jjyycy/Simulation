clc
clear
% a
r=0.05;
S0=100;
mu=0.1; % no use
sig=0.1;
T=1;
N=52;
K=100;
n=1000;
Sastore=zeros(N+1,n);
Sastore(1,:)=S0*ones(1,n); %initial value
delta=T/N;

for j=2:(N+1) 
    Sastore(j,:)=Sastore(j-1,:).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*randn(1,n));
end

Ca=zeros(1,n);
for i=1:n
    Ca(i)=max(1/N*sum(Sastore(2:(N+1),i))-K,0)*exp(-r*T);
end

Cbara=mean(Ca)
stderra=std(Ca)/sqrt(n)

% b
Sbstore1=zeros(N+1,n);
Sbstore1(1,:)=S0*ones(1,n); %initial value
Sbstore2=zeros(N+1,n);
Sbstore2(1,:)=S0*ones(1,n); %initial value
for j=2:(N+1)
    z=randn(1,n);
    Sbstore1(j,:)=Sbstore1(j-1,:).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*z);
    Sbstore2(j,:)=Sbstore2(j-1,:).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*(-z));
end

Cb=zeros(1,n);
for i=1:n
    Cb(i)=1/2*(max(1/N*sum(Sbstore1(2:(N+1),i))-K,0)+max(1/N*sum(Sbstore2(2:(N+1),i))-K,0))*exp(-r*T);
end
Cbarb=mean(Cb)
stderrb=std(Cb)/sqrt(n)

% c
Y=Ca;
X=Sastore(N+1,:);
a=-corr(X',Y')*std(Y)/std(X);
Cbarc=Cbara+a*(mean(X)-S0*exp(r*T))
stderrc=std(Y)/sqrt(n)*sqrt(1-(corr(X',Y'))^2)

% d
Yd=Ca;
Xd=zeros(1,n);
for i=1:n
    Xd(i)=(max((prod(Sastore(2:(N+1),i)))^(1/N)-K,0))*exp(-r*T);
end

% check that rho of Xd and Yd is huge
% rho=corr(Xd',Yd')
ad=-corr(Xd',Yd')*std(Yd)/std(Xd);

% calculate real EXd using approximation:
% sig=sig/sqrt(s)
% d=r/2+sig^2/12
EXd=european_call_div(S0, K, r, sig/sqrt(3), T, r/2+sig^2/12);

Cbard=Cbara+ad*(mean(Xd)-EXd)
stderrd=std(Yd)/sqrt(n)*sqrt(1-(corr(Xd',Yd'))^2)

% e
% calculate real EXe using:
sig_e=sig*sqrt((N+1)*(2*N+1)/(6*N^2));
d_e=r*(N-1)/(2*N)+sig^2*((N+1)*(N-1)/(12*N^2));
EXe=european_call_div(S0, K, r, sig_e, T, d_e);
Cbare=Cbara+ad*(mean(Xd)-EXe)
stderre=std(Yd)/sqrt(n)*sqrt(1-(corr(Xd',Yd'))^2)
