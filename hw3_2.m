clc
clear
% a
r=0.05;
S0=100;
mu=0.1;
sig=0.2;
T=1;
K=160; % or 140, 160
n=10000;

Z=randn(1,n);
STa=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Z);
Ca=exp(-r*T)*max(STa-K,0);
Cabar=mean(Ca)
stderra=std(Ca)/sqrt(n)

% b
Pb=exp(-r*T)*max(K-STa,0);
% Put-Call Parity:
% C-P=S0-K*exp(-r*T)
Cb=Pb+S0-K*exp(-r*T);
Cbbar=mean(Cb)
stderrb=std(Cb)/sqrt(n)

% c
Y=Cb;
X=STa;
a=-corr(X',Y')*std(Y)/std(X);
Ccbar=Cbbar+a*(mean(X)-S0*exp(r*T))
stderrc=std(Y)/sqrt(n)*sqrt(1-(corr(X',Y'))^2)

% d
L=(log(K/S0)-(r-0.5*sig^2)*T)/(sig*sqrt(T));
Ud=rand(1,n);
Xd=norminv(Ud.*(1-normcdf(L))+normcdf(L),0,1);
STd=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Xd);
Cd=exp(-r*T)*(STd-K)*(1-normcdf(L));
Cdbar=mean(Cd)
stderrd=std(Cd)/sqrt(n)

