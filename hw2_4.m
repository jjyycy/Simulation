n=30; % or 3, 5, 10, 30
m=1000;
x=zeros(1,m);
for i=1:m
    y=tan(rand(1).*pi-pi./2); % cauchy dist.
    u=rand(1);
    % C=sqrt(2*pi/e)
    while (u>1/sqrt(2*pi/exp(1))*gamma((1+n)/2)/sqrt((n*pi)*gamma(n/2))*(1+y.^2./n).^(n+1)/2/(1/(pi*(1+y.^2))))
        y=tan(rand(1).*pi-pi./2);
        u=rand(1);
    end
    x(i)=y;
end
qqplot(x);
%xlim([-2 2]);