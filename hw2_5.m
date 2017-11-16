% a
n=1000;
Z=randn(1,n);
S0=100;
T=1;
r=0.05;
sig=0.1;
K=105; % or 100, 105
S=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Z);
C=max(S-K,zeros(1,n))*exp(-r*T);
Cbar=mean(C)
stderr=sqrt(sum((C-Cbar).^2)/1000/999)


% real value
dplus=(log(S0/K)+(r+sig^2/2)*T)/(sig*sqrt(T));
dminus=dplus-sig*sqrt(T);
realC=S0*normcdf(dplus)-normcdf(dminus)*K*exp(-r*T)



% b
Z2=-Z;
S2=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Z2);
C2=max(S2-K,zeros(1,n))*exp(-r*T);
newC=(C+C2)./2;
newCbar=mean(newC)
newstderr=sqrt(sum((newC-newCbar).^2)/1000/999)


