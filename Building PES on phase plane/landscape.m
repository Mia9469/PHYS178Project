close all;clear
%由sim_trial.py生成
load('hist_s30_0.mat')

dist = 2;
s1 = S(1,:);
s2 = S(2,:);


% 坐标范围
X = -0.5:0.01:1;
Y = -0.5:0.01:1;
l = length(X); 
Z = zeros( l, l );

if dist==2
    p1=mean(s1);
    s1a = s1(min(s1>p1,s2<p1));
    s2a = s2(min(s1>p1,s2<p1));
    s1b = s1(min(s1<p1,s2>p1));
    s2b = s2(min(s1<p1,s2>p1));
    sigma = [std(s2a) 0; 0 std(s1a)];
    mu = [mean(s1a) mean(s2a)];
    p0=mean(mu);
    sigmab = [std(s2b) 0; 0 std(s1b)];
    mub = [mean(s1b) mean(s2b)];
    
    Zb = Z;
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma^(-2) * (p - mu)')/2; 
            Zb( row, col ) = - ((p - mub) * sigmab^(-2) * (p - mub)')/2;
        end
    end
 
     Z = exp(Z) / ( 2 * pi * det(sigma));
     Zb = exp(Zb) / ( 2 * pi * det(sigmab));
    Z = (Z*length(s1a) + Zb*length(s1b))/length(s1)*1e-4;
elseif dist == 3
    p1=max(s1)/2;
    p2=max(s2)/2;
    s1a = s1(min(s1<p1,s2>p2));
    s2a = s2(min(s1<p1,s2>p2));
    s1b = s1(min(s1<p1,s2<p2));
    s2b = s2(min(s1<p1,s2<p2));
    s1c = s1(min(s1>p1,s2<p2));
    s2c = s2(min(s1>p1,s2<p2));
    sigma = [std(s1a) 0; 0 std(s2a)];
    mu = [mean(s1a) mean(s2a)];
    sigmab = [std(s1b) 0; 0 std(s2b)];
    mub = [mean(s1b) mean(s2b)];
    sigmac = [std(s1c) 0; 0 std(s2c)];
    muc = [mean(s1c) mean(s2c)];
    
    Zb = Z;
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma^(-2) * (p - mu)')/2; 
            Zb( row, col ) = - ((p - mub) * sigmab^(-2) * (p - mub)')/2;
            Zc( row, col ) = - ((p - muc) * sigmac^(-2) * (p - muc)')/2;
        end
    end
 
     Z = exp(Z) / ( 2 * pi * det(sigma));
     Zb = exp(Zb) / ( 2 * pi * det(sigmab));
     Zc = exp(Zc) / ( 2 * pi * det(sigmac));
     Z =  (Z*length(s1a) + Zb*length(s1b)+ Zc*length(s1c))/length(s1)*1e-4;

else
    sigma =  [std(s1) 0; 0 std(s2)];
%     sigma = diag([ 0.018 0.065]);
%     sigma = diag([ 0.036 0.13]);
    mu = [mean(s1) mean(s2)];
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma^(-2) * (p - mu)')/2; 
        end
    end
 
    Z = exp(Z) / ( 2 * pi * det(sigma))*1e-4;
end

P = -log(Z);
 
% 三维曲面
figure, surf(X,Y, Z');
xlabel("S2");ylabel("S1");zlabel("Possibility");;xlim([0,1]);ylim([0,1])
figure, surf(X,Y, P');
xlabel("S2");ylabel("S1");zlabel("Potential");
xlim([-0.2,1]);ylim([-0.2,1])
zlim([0,300])
shading interp

figure, surf(edges(2,2:end), edges(1,2:end), hist);
xlabel("S2");ylabel("S1");zlabel("Possibility");

% figure, histogram2(s1,s2,100,'DisplayStyle','tile');
% xlabel("S2");ylabel("S1");
sum(sum(Z))*0.01^2


% for i=1:4
%     test(i,:) =  [std(s1(1:100000+i*100000)) std(s2(1:100000+i*100000))];
% end