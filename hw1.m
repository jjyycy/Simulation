Q = [-.1154, .1019, .0083, .0020, .0031, 0, 0, 0;
    .0091, -.1043, .0787, .0105, .0030, .0030, 0, 0;
    .0010, .0309, -.1172, .0688, .0107, .0048, 0, .0010;
    .0007, .0047, .0713, -.1711, .0701, .0174, .0020, .0049;
    .0005, .0025, .0089, .0813, -.2530, .1181, .0144, .0273;
    0, .0021, .0034, .0073, .0568, -.1928, .0479, .0753;
    0,0,.0142, .0142, .0250, .0928, -.4318, .2856;
    0, 0, 0, 0, 0, 0, 0, 0];
M = [0,.8838, .0720, .0173, .0269, 0, 0, 0;
    .0872, 0, .7545, .1007, .0288, .0288, 0,0;
    .0085, .2637, 0, .5870, .0913, .0410, 0, .0085;
    .0041, .0275, .4167, 0, .4097, .1017, .0117, .0286;
    .0020, .0099, .0352, .3213, 0, .4668, .0569, .1079;
    0, .0109, .0176, .0379, .2946, 0, .2484, .3906;
    0, 0, .0329, .0329, .0579, .2149, 0, .6614;
    0,0,0,0,0,0,0,1];
lam = [.1154, .1043, .1172, .1711, .2530, .1929, .4318, .0001];

% (a)
% S1. create aliasing table
AliasTable=zeros(56,3);
% aliasing table for state k : row 1+(k-1)*7 ~ 7+(k-1)*7
for k=1:8
    one2seven=[1;2;3;4;5;6;7;8];
    L=[one2seven,transpose(M(k,:))];
    N=size(L,1);
    L(:,2)=(N-1)*L(:,2);
    T=zeros(N-1,3);
    for i =1:N-1
        L=sortrows(L(1:N-i+1,:),2);
        T(i,1)=L(1,2);
        T(i,2)=L(1,1);
        T(i,3)=L(N-i+1,1);
        L(N-i+1,2)=L(N-i+1,2)-1+L(1,2);
        for j=2:N-i+1
            L(j-1,1)=L(j,1);
            L(j-1,2)=L(j,2);
        end
    end
    AliasTable((1+(k-1)*7):(7+(k-1)*7),:)=T;
end

% S2. loop through
n=1000000; % generate 1000,000 paths

count=zeros(8,8);

for k=1:8 % states
    for i=1:n
        % before time expires
        %clock=(-1/lam(k)*log(rand));
        % (b)
        clock=1/2*((-1/lam(k)*log(rand))+(-1/lam(k)*log(rand)));
        state=k;
        while clock<5
            % random select next state
            T=AliasTable((1+(state-1)*7):(7+(state-1)*7),:);
            P=transpose(T(:,1));
            X=transpose(T(:,2));
            A=transpose(T(:,3));
            V=(8-1)*rand(1);
            I=ceil(V);
            W=I-V;
            Y=(W<=P(I));
            X=X(I)*Y+A(I)*(1-Y);
            state=X;
            %clock=clock+(-1/lam(state)*log(rand));  % holding time
            % (b)
            clock=clock+1/2*((-1/lam(k)*log(rand))+(-1/lam(k)*log(rand)));
        end
        count(state,k)=count(state,k)+1;
    end 
end

transmat=transpose(count)/n;
%std. error matrix
stderr=sqrt(n*transmat.*(1-transmat))/n;

%P5_2=expm(Q*5)



