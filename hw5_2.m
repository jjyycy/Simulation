clc
clear

% % a
% T=0.25;
% N=30;
% n=1000;
% S0=50;
% mu=0.15;
% sig=0.25;
% r=0.1;
% K=S0;
% 
% % standard Monte Carlo
% delta=T/N;
% S_a1=zeros(N+1,n);
% S_a1(1,:)=S0*ones(1,n);
% 
% for j=2:(N+1)
%     z=randn(1,n);
%     S_a1(j,:)=S_a1(j-1,:)+r*delta*S_a1(j-1,:)+sqrt(delta)*sig*S_a1(j-1,:).*z;
% end
% 
% C_a1=exp(-r*T)*max(max(S_a1)-S0,0);
% Cbar_a1=mean(C_a1)
% Cstderr_a1=std(C_a1)/sqrt(n)
% 
% M=zeros(N,n);
% % brownian bridge
% for j=1:N
%     % consider time period j*delta,(j+1)*delta
%     b=(S_a1(j+1,:)-S_a1(j,:))./(sig.*S_a1(j,:)); % B_end
%     u=rand(1,n);
%     maxB=(b+sqrt(b.^2-2*delta*log(1-u)))/2;
%     M(j,:)=S_a1(j,:)+sig*S_a1(j,:).*maxB;
% end
% 
% C_a2=exp(-r*T)*max(max(M)-S0,0);
% Cbar_a2=mean(C_a2)
% Cstderr_a2=std(C_a2)/sqrt(n)
% 
% % b
% %comment all above
% T=0.25;
% N=30;
% n=1000;
% S0=50;
% mu=0.15;
% sig=0.25;
% r=0.1;
% 
% % standard Monte Carlo
% delta=T/N;
% S_b1=zeros(N+1,n);
% S_b1(1,:)=S0*ones(1,n);
% 
% for j=2:(N+1)
%     z=randn(1,n);
%     S_b1(j,:)=S_b1(j-1,:)+r*delta*S_b1(j-1,:)+sqrt(delta)*sig*S_b1(j-1,:).*z;
% end
% 
% C_b1=exp(-r*T)*max(max(S_b1)-S_b1(N+1,:),0);
% Cbar_b1=mean(C_b1)
% Cstderr_b1=std(C_b1)/sqrt(n)
% 
% M=zeros(N,n);
% % brownian bridge
% for j=1:N
%     % consider time period j*delta,(j+1)*delta
%     b=(S_b1(j+1,:)-S_b1(j,:))./(sig.*S_b1(j,:)); % B_end
%     u=rand(1,n);
%     maxB=(b+sqrt(b.^2-2*delta*log(1-u)))/2;
%     M(j,:)=S_b1(j,:)+sig*S_b1(j,:).*maxB;
% end
% 
% C_b2=exp(-r*T)*max(max(M)-S_b1(N+1,:),0);
% Cbar_b2=mean(C_b2)
% Cstderr_b2=std(C_b2)/sqrt(n)

% c
% comment all above
T=0.25;
N=30;
n=100000;
S0=50;
mu=0.15;
sig=0.5;
r=0.1;
K=50;
H=45;

% standard Monte Carlo
delta=T/N;
S_c1=zeros(N+1,n);
S_c1(1,:)=S0*ones(1,n);

for j=2:(N+1)
    z=randn(1,n);
    S_c1(j,:)=S_c1(j-1,:)+r*delta*S_c1(j-1,:)+sqrt(delta)*sig*S_c1(j-1,:).*z;
end

C_c1=zeros(1,n);
minS=min(S_c1);

for i=1:n
    if (minS(i)>H)
        C_c1(i)=exp(-r*T)*max(S_c1(N+1,i)-K,0);
    end
end
Cbar_c1=mean(C_c1)
Cstderr_c1=std(C_c1)/sqrt(n)

Mm=zeros(N,n);
% brownian bridge
for j=1:N
    % consider time period j*delta,(j+1)*delta
    b=(S_c1(j+1,:)-S_c1(j,:))./(sig.*S_c1(j,:)); % B_end
    u=rand(1,n);
    minB=(b-sqrt(b.^2-2*delta*log(1-u)))/2;
    Mm(j,:)=S_c1(j,:)+sig*S_c1(j,:).*minB;
end

C_c2=zeros(1,n);
minMm=min(Mm);
for i=1:n
    if (minMm(i)>H)
        C_c2(i)=exp(-r*T)*max(S_c1(N+1,i)-K,0);
    end
end
Cbar_c2=mean(C_c2)
Cstderr_c2=std(C_c2)/sqrt(n)

