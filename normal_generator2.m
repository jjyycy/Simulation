% a
seed=1;
rand('seed',seed);
n=100;
X = zeros(1,n);
for i=1:n
    s=sign(rand(1)-0.5);
    u1=rand(1);
    u2=rand(1);
    x1=-log(u1);
    while (u2>exp(-0.5*(x1-1)^2))
        u1=rand(1);
        u2=rand(1);
        x1=-log(u1);
    end;
    X(i)=x1*s;
end;

qqplot(X);

% b
Y = zeros(1,n);
for i=1:n
    x=rand(1);
    Y(i)=0+(x^0.1349-(1-x)^0.1349)/0.1975;
end

qqplot(Y);

% c
Z= zeros(1,n);
for i=1:n
    u=rand(1);
    z=randn(1);
    if (u<=0.82)
        Z(i)=0.6*z;
    else
        Z(i)=1.98*z;
    end      
end
    
qqplot(Z);


