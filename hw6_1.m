clc
clear
S0=40;
K=40; % in this case will not exercise at time 0
sig=0.2;
r=0.06;
T=2; % or 2
N=50;
n=50000;
delta=T/N;

% First generate stock price paths
S=zeros(n,N);
San=zeros(n,N); %store antithetic values of S
S0v=ones(n,1)*S0;

for j=1:N
    z=randn(n,1);
    if j==1   
        S(:,j)=S0v.*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*z);
        San(:,j)=S0v.*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*(-z));
    else
        S(:,j)=S(:,j-1).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*z);
        San(:,j)=San(:,j-1).*exp((r-1/2*sig^2)*delta+sig*sqrt(delta).*(-z));
    end
end
S=[S;San];% first 50000 rows: original; second 50000 rows: antithetic
P = max(K-S(:,50),0); % Payoff at final time

% Cash flow table
for i = (N-1):-1:1
    Si = S(:,i);
    itmP = find(K-Si>0); % ITM path label
    otmP = find(K-Si<=0);
    X = Si;
    Y = P*exp(-r*delta); % Discounted payoffs (a reduced-dim column)
    Y(otmP) = 0;
    XMat = [ones(size(X)),Laguerre(X,0),Laguerre(X,1),Laguerre(X,2)]; % MODIFIED
    XMat(otmP,:)=0;
    beta = XMat\Y;
    C = XMat*beta; % Value of continuation
    E = K-X; % Value of immediate exercise
    exP = find(C<E); % immediate exercise
    P(exP) = E(exP);
    rest = setdiff(1:n,exP);
    P(rest) = P(rest)*exp(-r*delta); % continue
end
P1=P(1:n)*exp(-r*delta);
Pa=P((n+1):2*n)*exp(-r*delta);
newP=(P1+Pa)/2;
P_bar = mean(newP)
P_stderr = std(newP)/sqrt(n)

