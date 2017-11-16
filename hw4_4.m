clc
clear
S0=95;
K=100;
r=0.05;
sig=0.2;
T=1;
n=10000;
h=0.0001;

% a
d_minus=1/(sig*sqrt(T))*(log(S0/K)+(r-1/2*sig^2)*T);
delta_a=exp(-r*T)*normpdf(d_minus)*1/(sig*sqrt(T))*1/S0

% b
z=randn(n,1);
S=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*z);
C=exp(-r*T)*(S>K); % pays 1 or 0
Sh=(S0+h)*exp((r-1/2*sig^2)*T+sig*sqrt(T).*z);
Ch=exp(-r*T)*(Sh>K);
delta_b=(Ch-C)./h;
delta_b_bar=mean(delta_b)
delta_b_stderr=std(delta_b)/sqrt(n)

% c
z=randn(n,1);
S=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*z);
delta_c=exp(-r*T)*(S>K)*1/(S0*sig^2*T).*(log(S./S0)-(r-0.5*sig^2)*T);
delta_c_bar=mean(delta_c)
delta_c_stderr=std(delta_c)/sqrt(n)