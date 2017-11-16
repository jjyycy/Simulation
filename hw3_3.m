clear
clc

% a
r=0.05;
S0=95;
sig=0.15;
T=0.25;
m=50;
% (H,K)=(94,96) or (90,96) or (85,96) or (90, 106)
H=90;
K=106;
n=100000;
Sastore=zeros(m+1,n);
Sastore(1,:)=S0*ones(1,n); %initial value
delta=T/m;

for j=2:(m+1) 
    Sastore(j,:)=Sastore(j-1,:).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*randn(1,n));
end

Ca=zeros(1,n);
for i=1:n
    if (min(Sastore(:,i))<H)
        if (Sastore(m+1,i)-K>0)
            Ca(i)=1*exp(-r*T);
        else
            Ca(i)=0;
        end
    else
        Ca(i)=0;
    end
end
Cbara=mean(Ca)
stderra=std(Ca)/sqrt(n)

% b
% (H,K)=(94,96) or (90,96) or (85,96) or (90, 106)
Sbstore=zeros(m+1,n);
Sbstore(1,:)=S0*ones(1,n); %initial value
Cb=zeros(1,n);

for i=1:n
    for j=2:(m+1)
        Sastore(j,i)=Sastore(j-1,i).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*randn(1));
        if (Sastore(j,i)<H)
            Cb(i)=exp(-delta*(j-1))*digital_call(Sastore(j,i), K, r, sig, T-(j-1)*delta);
            break
        end
    end
end
Cbarb=mean(Cb)
Cbarb*10000
stderrb=std(Cb)/sqrt(n)

(stderra/stderrb)^2

