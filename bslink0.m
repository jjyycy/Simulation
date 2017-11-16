function Y=bslink0(m,X)
S0=100;
r=0.05;
sig=0.20;
T=1;
K=160; % or 140, 160
G=exp(-r*T)*max(S0*exp((r-sig^2/2)*T+X*sig*sqrt(T))-K,0);
W=exp(-m*X+m^2/2);
%Y=G.*sqrt(W);
Y=((exp(-r*T)*max(S0*exp((r-sig^2/2)*T+X*sig*sqrt(T))-K,0)).*sqrt(exp(-m*X+m^2/2)));