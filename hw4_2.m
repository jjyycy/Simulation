clc
clear
% a
alpha=0.2;
sig=0.1;
b=0.05;
r0=0.04;
n=1000;
N=50;
T=1;
delta=T/N;

r=zeros(n,N+1);
r(:,1)=r0*ones(n,1);
sumr=zeros(n,1);
for i=2:(N+1)
    d=4*alpha*b/(sig^2);
    lam=4*alpha*exp(-alpha*delta)/(sig^2*(1-exp(-alpha*delta))).*r(:,i-1);
    chi=random('ncx2',d,lam);
    r(:,i)=(sig^2*(1-exp(-alpha*delta)))/(4*alpha).*chi; 
    sumr=sumr+r(:,i);
end

ZCB=exp(-delta*sumr)
ZCB_bar=mean(ZCB)
ZCB_stderr=std(ZCB)/sqrt(n)

% b
t=1;
del=1/12;
L=1;
R=0.05;

payoff=L*del*max(0,r(:,N+1)-R);
caplet=exp(-delta.*sumr).*payoff;
caplet_bar=mean(caplet)
caplet_stderr=std(caplet)/sqrt(n)