clc
clear
S0=100;
r=0.05;
sig=0.20;
T=1;
K=120; % or 140, 160
n=10000;

% a
% B-S formula
dplus=1/(sig*sqrt(T))*(log(S0/K)+(r+1/2*sig^2)*T);
dminus=1/(sig*sqrt(T))*(log(S0/K)+(r-1/2*sig^2)*T);
BS_price = S0*normcdf(dplus)-normcdf(dminus)*K*exp(-r*T)

% standard MC
S_stdMC=S0*ones(n,1).*exp((r-0.5*sig^2)*T+sig*sqrt(T).*randn(n,1));
C_stdMC=exp(-r*T).*max(S_stdMC-K,0);
C_stdMC_bar=mean(C_stdMC)
C_stdMC_stderr=std(C_stdMC)/sqrt(n)

% GHS
% solve for m I:
syms x positive
S_T_m=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T)*x);
eqn = x == S_T_m*sig*sqrt(T)/(S_T_m-K);
xsolve = vpasolve(eqn,3)
m=double(xsolve)

% solve for m II:
solve_fun = @(x) sig*sqrt(T)/(1-K/(S0*exp((r-1/2*sig^2)*T+sig*sqrt(T)*x)))-x;
low = 1.0;
high = 2.0;
tolerance = .00001;
x = bisection(solve_fun, low, high, tolerance); 




% importance sampling
z=randn(n,1);
S_IS=S0*ones(n,1).*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(z+m));
C_IS=exp(-r*T).*max(S_IS-K,0).*exp(-m*z-0.5*m^2);
C_IS_bar=mean(C_IS)
C_IS_stderr=std(C_IS)/sqrt(n)

% b
K=160; % or 140, 160
M=100000;
ind=randn(M,1);
dep=zeros(M,1);
mtrial=1;
%bslink = @(m, X)((exp(-r*T)*max(S0*exp((r-sig^2/2)*T+X*sig*sqrt(T))-K,0)).*sqrt(exp(-m*X+m^2/2)));
mhat=nlinfit(ind,dep,@bslink0,mtrial)
z=randn(n,1);
S_C=S0*ones(n,1).*exp((r-0.5*sig^2)*T+sig*sqrt(T)*(z+mhat));
C_C=exp(-r*T).*max(S_C-K,0).*exp(-mhat*z-0.5*mhat^2);
C_C_bar=mean(C_C);
C_C_stderr=std(C_C)/sqrt(n)

% c
K=160; % or 140, 160
bslink = @(m, X)((exp(-r*T)*max(S0*exp((r-sig^2/2)*T+X*sig*sqrt(T))-K,0)).*sqrt(exp(m(2))*exp(-0.5*(X.^2-((X-m(1))/exp(m(2))).^2))));
ind = randn(M,1);
dep = zeros(M,1);
mtrial = [1,1];
mhat = nlinfit(ind,dep,bslink,mtrial);
m = mhat(1)
s = exp(mhat(2))
