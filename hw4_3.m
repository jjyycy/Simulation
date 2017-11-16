clc
clear
r=0.1;
K=100;
q=0.03;
sig=0.25;
T=0.2;
n=10000;
h=0.0001;
S0=90;
%European call with dividend

%% resimulation
% without control
z=randn(n,1);
S=S0*exp((r-q-1/2*sig^2)*T+sig*sqrt(T).*z);
C=exp(-r*T)*max(S-K,0);
Sh=(S0+h)*exp((r-q-1/2*sig^2)*T+sig*sqrt(T).*z);
Ch=exp(-r*T)*max(Sh-K,0);
delta=(Ch-C)./h;
delta_bar=mean(delta)
delta_stderr=std(delta)/sqrt(n)

Sv=S0*exp((r-q-1/2*(sig+h)^2)*T+(sig+h)*sqrt(T).*z);
Cv=exp(-r*T)*max(Sv-K,0);
vega=(Cv-C)./h;
delta_bar=mean(vega)
delta_stderr=std(vega)/sqrt(n)

% with control
Yd=delta;
X=S;
EX=exp((r-q)*T)*S0*ones(n,1);
deltac=Yd+(-corr(Yd,X)*std(Yd)/std(X))*(mean(X)-EX);
deltac_bar=mean(deltac)
deltac_stderr=std(deltac)/sqrt(n)*sqrt(1-(corr(Yd,X))^2)

Yv=vega;
vegac=Yv+(-corr(Yv,X)*std(Yv)/std(X))*(mean(X)-EX);
vegac_bar=mean(vegac)
vegac_stderr=std(vegac)/sqrt(n)*sqrt(1-(corr(Yv,X))^2)

%% pathwise
z=randn(n,1);
S=S0*exp((r-q-1/2*sig^2)*T+sig*sqrt(T).*z);

% without control
delta1=zeros(n,1);
vega1=zeros(n,1);
for i=1:n
    if (S(i)>=K)
        delta1(i)=exp(-r*T)*S(i)/S0;
        vega1(i)=exp(-r*T)*S(i)/sig*(log(S(i)/S0)-(r-q+1/2*sig^2)*T);
    end
end

delta1_bar=mean(delta1)
delta1_stderr=std(delta1)/sqrt(n)

vega1_bar=mean(vega1)
vega1_stderr=std(vega1)/sqrt(n)

% with control
Y1d=delta1;
X=S;
EX=exp((r-q)*T)*S0*ones(n,1);
add=-corr(Y1d,X)*std(Y1d)/std(X);
delta1c=Y1d+add*(mean(X)-EX);
delta1c_bar=mean(delta1c)
delta1c_stderr=std(delta1c)/sqrt(n)*sqrt(1-(corr(Y1d,X))^2)

Y1v=vega1;
adv=-corr(Y1v,X)*std(Y1v)/std(X);
vega1c=Y1v+adv*(mean(X)-EX);
vega1c_bar=mean(vega1c)
vega1c_stderr=std(vega1c)/sqrt(n)*sqrt(1-(corr(Y1v,X))^2)

%% loglikelihood
% without control
z=randn(n,1);
S=S0*exp((r-q-1/2*sig^2)*T+sig*sqrt(T).*z);
delta2=exp(-r*T)*max(S-K,0)*1/(S0*sig^2*T).*(log(S./S0)-(r-q-0.5*sig^2)*T);
delta2_bar=mean(delta2)
delta2_stderr=std(delta2)/sqrt(n)

d=(log(S./S0)-(r-q-sig^2/2)*T)/(sig*sqrt(T));
dd_over_dsig=(log(S0./S)+(r-q+sig^2)*T)/(sig^2*sqrt(T));
dlogg_over_dsig=-d.*dd_over_dsig-1/sig;
vega2=exp(-r*T).*max(S-K,0).*dlogg_over_dsig;
vega2_bar=mean(vega2)
vega2_stderr=std(vega2)/sqrt(n)

% with control
Y2d=delta2;
X=S;
EX=exp((r-q)*T)*S0*ones(n,1);
ad2=-corr(Y2d,X)*std(Y2d)/std(X);
delta2c=Y2d+ad2*(mean(X)-EX);
delta2c_bar=mean(delta2c)
delta2c_stderr=std(delta2c)/sqrt(n)*sqrt(1-(corr(Y2d,X))^2)

Y2v=vega2;
adv2=-corr(Y2v,X)*std(Y2v)/std(X);
vega2c=Y2v+adv*(mean(X)-EX);
vega2c_bar=mean(vega2c)
vega2c_stderr=std(vega2c)/sqrt(n)*sqrt(1-(corr(Y2v,X))^2)

