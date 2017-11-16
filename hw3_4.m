clc
% a
r=0.10;
S0=100;
sig=0.30;
T=0.2;
N=50; % or 50
H=95;
K=100;
n=100000;
Sastore=zeros(N+1,n);
Sastore(1,:)=S0*ones(1,n); %initial value
delta=T/N;

for j=2:(N+1) 
    Sastore(j,:)=Sastore(j-1,:).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*randn(1,n));
end

Ca=zeros(1,n);
for i=1:n
    if (min(Sastore(:,i))<H)
        Ca(i)=exp(-r*T)*max(Sastore(j,i)-K,0);
    end
end
Cbara=mean(Ca)
stderra=std(Ca)/sqrt(n)

% b
Sbstore=zeros(N+1,n);
Sbstore(1,:)=S0*ones(1,n); %initial value
Cb=zeros(1,n);

for i=1:n
    for j=2:(N+1)
        Sbstore(j,i)=Sbstore(j-1,i).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*randn(1));
        if (Sbstore(j,i)<H)
            ST=Sbstore(j,i)*exp((r-0.5*sig^2)*(T-(j-1)*delta)+sig*sqrt(T-(j-1)*delta)*randn(1));
            Cb(i)=exp(-r*T)*max(ST-K,0);
            break
        end
    end
end
Cbarb=mean(Cb)
stderrb=std(Cb)/sqrt(n)

% continuous time
lam=(r+sig^2/2)/(sig^2);
y=(log(H^2/(S0*K)))/(sig*sqrt(T))+lam*sig*sqrt(T);
Cdi=S0*(H/S0)^(2*lam)*normcdf(y)-K*exp(-r*T)*(H/S0)^(2*lam-2)*normcdf(y-sig*sqrt(T))

% c
N=50; % or 50
theta=-0.30; % or -0.30
Scstore=zeros(N+1,n);
Scstore(1,:)=S0*ones(1,n); %initial value
Cc=zeros(1,n);

for i=1:n
    multiplier=1;
    for j=2:(N+1)
        z=randn(1);
        Scstore(j,i)=Scstore(j-1,i).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*(z+theta)); % (z+theta) ~ N(theta,1)
        multiplier=multiplier*exp(-theta*z-0.5*theta^2);
        if (Scstore(j,i)<H)
            Cc(i)=exp(-delta*(j-1))*european_call_div(Scstore(j,i), K, r, sig, T-(j-1)*delta, 0)*multiplier;
            break
        end
    end
end
Cbarc=mean(Cc)
stderrc=std(Cc)/sqrt(n)