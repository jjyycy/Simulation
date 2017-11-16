clc
clear

S10=100;
S20=100;
K=100;
r=0.1;
sig1=0.3;
sig2=0.3;
rho=0.5;
T=0.2;
H=95;
n=10000;
N=50;


% a
S1store=zeros(N+1,n);
S2store=zeros(N+1,n);
S1store(1,:)=S10*ones(1,n); %initial value
S2store(1,:)=S20*ones(1,n); %initial value
delta=T/N;

for j=2:(N+1)
    z1=randn(1,n);
    z2=rho*z1+sqrt(1-rho^2)*randn(1,n);
    S1store(j,:)=S1store(j-1,:).*exp((r-1/2*sig1^2)*delta+sig1*sqrt(delta).*z1);
    S2store(j,:)=S2store(j-1,:).*exp((r-1/2*sig2^2)*delta+sig2*sqrt(delta).*z2);
end

minS2=min(S2store);
C1=zeros(1,n);
for i=1:n
    if (minS2(i)>95)
        C1(i)=exp(-r*T)*max(S1store(N+1,i)-K,0);
    end
end

C1_bar=mean(C1)
C1_stderr=std(C1)/sqrt(n)

% b
Mm=zeros(N,n);
for j=1:N
    % consider time period j*delta,(j+1)*delta
    b=(S2store(j+1,:)-S2store(j,:))./(sig2.*S2store(j,:)); % B_end
    u=rand(1,n);
    minB=(b-sqrt(b.^2-2*delta*log(1-u)))/2;
    Mm(j,:)=S2store(j,:)+sig2*S2store(j,:).*minB;
end

minMm=min(Mm);
C2=zeros(1,n);

for i=1:n
    if (minMm(i)>H)
        C2(i)=exp(-r*T)*max(S1store(N+1,i)-K,0);
    end
end
C2_bar=mean(C2)
C2_stderr=std(C2)/sqrt(n)