% clc
% clear
% N=5;
% T=5;
% r=0.04;
% s=0.01;
% R=0.35;
% lam=s/(1-R);
% corr=0;   % or 0.2, 0.4, 0.6, 0.8
% itd=5;  % or 2, 3, 4, 5
% n=100000;   
% 
% %Gaussian Copula
% Sigma=corr*ones(5)+(1-corr)*eye(5);
% A=chol(Sigma,'lower');
% z=randn(5,n);
% y=A*z; % or y=z if rho=1
% 
% U=normcdf(y);
% % exp(lam) cdf: 1-exp(-lam*x)
% X=-1/lam*log(1-U);
% 
% count=zeros(1,n);
% C=zeros(1,n);
% for i=1:n
%     count(i)=sum(X(:,i)<=5);
%     if count(i)>=itd
%         C(i)=exp(-r*T)*(1-R);
%     end
% end
% 
% C_bar=mean(C)
% C_stderr=std(C)/sqrt(n)


%% if rho=1
% comment all above
N=5;
T=5;
r=0.04;
s=0.01;
R=0.35;
lam=s/(1-R);
corr=1;
itd=5;  % or 2, 3, 4, 5
n=100000;   

%Gaussian Copula
z=randn(1,n);
y=[z;z;z;z;z];

U=normcdf(y);
% exp(lam) cdf: 1-exp(-lam*x)
X=-1/lam*log(1-U);

count=zeros(1,n);
C=zeros(1,n);
for i=1:n
    count(i)=sum(X(:,i)<=5);
    if count(i)>=itd
        C(i)=exp(-r*T)*(1-R);
    end
end

C_bar=mean(C)
C_stderr=std(C)/sqrt(n)

