%尝试使用未对角化的sigma矩阵
close all;clear
load('fixpoints.mat')
load('sigma_cor.mat')

fix=fix_30;
s=size(fix);
sig=sig_30;

dist=2;
if s(1)==5
    dist = 3;
    sig_c=sig_0c;
end
% 坐标范围
X = 0:0.01:1;
Y = 0:0.01:1;
l = length(X); 
Z = zeros( l, l );

if dist==2
    sigma = rot90(sig,2);
    mu = fix(1,:);
    sigmab = sig;
    mub = fix(3,:);
    
    Zb = Z;
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma.^(-2) * (p - mu)')/2; 
            Zb( row, col ) = - ((p - mub) * sigmab.^(-2) * (p - mub)')/2;
        end
    end
 
     Z = exp(Z) / ( 2 * pi * det(sigma));
     Zb = exp(Zb) / ( 2 * pi * det(sigmab));
    Z = (Z + Zb)/2*1e-4;
elseif dist == 3
    
    sigma =rot90(sig,2);
    mu = fix(1,:);
    sigmab = sig_c;
    mub = fix(3,:);
    sigmac = sig;
    muc = fix(5,:);
    
    Zb = Z;
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma.^(-2) * (p - mu)')/2; 
            Zb( row, col ) = - ((p - mub) * sigmab.^(-2) * (p - mub)')/2;
            Zc( row, col ) = - ((p - muc) * sigmac.^(-2) * (p - muc)')/2;
        end
    end
 
     Z = exp(Z) / ( 2 * pi * det(sigma));
     Zb = exp(Zb) / ( 2 * pi * det(sigmab));
     Zc = exp(Zc) / ( 2 * pi * det(sigmac));
     Z =  (Z + Zb+ Zc)/3*1e-4;
end

P = -log(Z);
 
% 三维曲面
figure, surf(X,Y, Z');
xlabel("S2");ylabel("S1");zlabel("Possibility");xlim([0,1]);ylim([0,1])
figure, surf(X,Y, P');
xlabel("S2");ylabel("S1");zlabel("Potential");
xlim([0.2,1]);ylim([0.2,1]);%zlim([-600,300])
title("Mu_0=30Hz")
shading interp
figure;contour(X,Y,P',100);
xlabel("S2");ylabel("S1");zlabel("Potential");

row=find(abs(X-fix(1,1))==min(abs(X-fix(1,1))));
col=find(abs(Y-fix(1,2))==min(abs(Y-fix(1,2))));
E=P(row,col)
E_w=min(min(P))
if s(1)==5
    c=find(abs(Y-fix(3,1))==min(abs(Y-fix(3,1))));
    s1=find(abs(X-fix(2,1))==min(abs(X-fix(2,1))));
    s2=find(abs(Y-fix(2,2))==min(abs(Y-fix(2,2))));
    [r,r2]=find(P==min(min(P(51:92,51:92))));
    E_s=P(s1,s2)
    E_c=P(c,c)
    E_lowc=P(r,r2)
%      save('E-12','E','E_c','E_lowc','r','E_w','E_s')
else
    s=find(abs(Y-fix(2,1))==min(abs(Y-fix(2,1))));
    E_s=P(s,s)
%       save('E-30','E','E_s','E_w')
end
sum(sum(Z))