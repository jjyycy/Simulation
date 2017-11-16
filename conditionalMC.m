clc
r=0.10;
S0=100;
sig=0.30;
T=0.2;
N=50; % or 50
H=95;
K=100;
n=100000;
delta=T/N;
C=zeros(1,n);
m=-0.3;

for i=1:n
    S=S0;
    t=0;
    RN=1;
    while (S>=H) && (t<T)
        z=randn(1);
        z_alt=z+m;
        t=t+delta;
        S=S*exp((r-1/2*sig^2)*delta+sig*sqrt(delta)*z_alt);
        RN=RN*exp(-m*z-1/2*m^2);
    end
    if t>=T
        C(i)=0;
    else
        C(i)=exp(-r*t)*european_call_div(S, K, r, sig, T-t, 0)*RN;
    end
end
C_bar=mean(C)
C_stderr=std(C)/sqrt(n)