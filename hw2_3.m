% % a
x=zeros(2,1000);
% rho=0;% or 0.4, 0.8
% x(1,:)=randn(1,1000);
% x(2,:)=rho.*x(1,:)+sqrt(1-rho^2).*randn(1,1000);
% scatter(x(1,:),x(2,:));
% 
% b
% y=zeros(2,1000);
% z=randn(2,1000); % Z with N(0,I) distribution
% rho=0.8; % or 0.4, 0.8
% y(1,:)=z(1,:);
% y(2,:)=rho*z(1,:)+sqrt(1-rho^2)*z(2,:);
% s=gamrnd(5/2,1/2,1,1000);
% T1=sqrt(5./s).*y(1,:);
% T2=sqrt(5./s).*y(2,:);
% scatter(T1,T2);

% % c
% rho=0.8; % or 0.4, 0.8
% z=randn(2,1000); % Z with N(0,1) distribution
% y(1,:)=z(1,:);
% y(2,:)=rho*z(1,:)+sqrt(1-rho^2)*z(2,:);
% U=normcdf(y);
% % exp(1) cdf: 1-exp(-x)
% X=-log(1-U);
% scatter(X(1,:),X(2,:));

% d
rho=0.8; % or 0.4, 0.8
z=randn(2,1000); % Z with N(0,1) distribution
y(1,:)=z(1,:);
y(2,:)=rho*z(1,:)+sqrt(1-rho^2)*z(2,:);
s=gamrnd(5/2,1/2,1,1000);
w=sqrt(5./([1;1]*s)).*y;
% t copula
U=tcdf(w,5);
% exp(1) cdf: 1-exp(-x)
X=-log(1-U);
scatter(X(1,:),X(2,:));
