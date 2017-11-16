% (a) normal generator
X100 = normrnd(0,1,[1,100]);
X1000 = normrnd(0,1,[1,1000]);
X10000 = normrnd(0,1,[1,10000]);
%qqplot(X100);
%qqplot(X1000);
%qqplot(X10000);

% (b) poor man's generator
Y100=zeros(1,100);
for i = 1:100
    sumu=0;
    for j=1:12
        sumu=sumu+rand;
    end
    Y100(i)=sumu-6;
end

Y1000=zeros(1,1000);
for i = 1:1000
    sumu=0;
    for j=1:12
        sumu=sumu+rand;
    end
    Y1000(i)=sumu-6;
end

Y10000=zeros(1,10000);
for i = 1:10000
    sumu=0;
    for j=1:12
        sumu=sumu+rand;
    end
    Y10000(i)=sumu-6;
end
%qqplot(Y100);
%qqplot(Y1000);
qqplot(Y10000);