clc
clear
S0=100;
sig=0.2;
T=1;
r=0.05;
K=100;
n=1000;

% a
U=rand(n,2);
Z=norminv(U);
S1=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Z(:,1));
S2=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*Z(:,2));
C_a=exp(-r*T)*max(1/2*(S1+S2)-K,0);
C_a_bar=mean(C_a)
C_a_stderr=std(C_a)/sqrt(n)

% b
% generate 10 points within each bin
C_b=[];
sqrerr_b=[];
for i=1:10
    for j=1:10
        % in bin (i-1)/10~i/10 , (j-1)/10~j/10)
        u=rand(10,2);
        u(:,1)=(i-1)/10+1/10.*u(:,1);
        u(:,2)=(j-1)/10+1/10.*u(:,2);
        z=norminv(u);
        s1=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*z(:,1));
        s2=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T).*z(:,2));
        C_temp=exp(-r*T)*max(1/2*(s1+s2)-K,0);
        C_b=[C_b;C_temp];
        % append square error for each bin
        sqrerr_b=[sqrerr_b;(std(C_temp))^2];
    end
end
C_b_bar=mean(C_b)
C_b_stderr=sqrt(1/100*sum(sqrerr_b))/sqrt(1000)

% c
C_c=[];
sqrerr_c=[];
for i=1:250
    % each bin: (i-1)/250~i/250
    nu=[1/sqrt(2);1/sqrt(2)];
    for k=1:4 % sample size n=4 in each bin
        u=(i-1)/250+1/250.*rand(1);
        x=norminv(u);
        Mu=x*nu;
        Sigma=eye(2)-nu*nu';
        z=mvnrnd(Mu,Sigma);
        s1=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T)*z(1));
        s2=S0*exp((r-1/2*sig^2)*T+sig*sqrt(T)*z(2));
        C_temp=exp(-r*T)*max(1/2*(s1+s2)-K,0);
        C_c=[C_c;C_temp];
    end
    % append square error for each bin
    sqrerr_c=[sqrerr_c;(std(C_c((i-1)*4+1:(i-1)*4+4)))^2];
end
C_c_bar=mean(C_c)
C_c_stderr=sqrt(1/250*sum(sqrerr_c))/sqrt(1000)