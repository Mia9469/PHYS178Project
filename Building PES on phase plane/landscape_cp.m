close all;clear
load('fixpoints_cp.mat')
load('sigma_cp.mat')

fix=fix_0;
s=size(fix);
sig=sig_0;
p=[0.5000, 0.5350, 0.5432, 0.5718, 0.5871, 0.6011, 0.6188, 0.6336, 0.6395];
por=p(2);
% ×ø±ê·¶Î§
X = -0.5:0.01:1;
Y = -0.5:0.01:1;
l = length(X); 
Z = zeros( l, l );
Zb = Z;

sigma = diag(sig(1,:));
mu = fix(1,:);
sigmab = diag(sig(2,:));
mub = fix(3,:);
    
    for row = 1 : 1 : l
        for col = 1 : 1 : l
            p=[X(row) Y(col)];
            Z( row, col ) = - ((p - mu) * sigma^(-2) * (p - mu)')/2; 
            Zb( row, col ) = - ((p - mub) * sigmab^(-2) * (p - mub)')/2;
        end
    end
 
     Z = exp(Z) / ( 2 * pi * det(sigma));
     Zb = exp(Zb) / ( 2 * pi * det(sigmab));
    Z = (Z.*por + Zb.*(1-por))*1e-4;
    P = -log(Z);
    
 
% ÈýÎ¬ÇúÃæ
figure, surf(X,Y, Z);
xlabel("S2");ylabel("S1");zlabel("Possibility");
figure, surf(X,Y, max(P,10^-100));hold on;
plot3(0,0,P(51,51)+20,'y*','LineWidth',5);
xlabel("S2");ylabel("S1");zlabel("Potential");
xlim([-0.2,1]);ylim([-0.2,1]);zlim([0,300])
title("Cprime=0%")
shading interp
figure;contour(X,Y,P,100);
xlabel("S2");ylabel("S1");zlabel("Potential");

row=find(abs(X-fix(1,1))==min(abs(X-fix(1,1))));
col=find(abs(Y-fix(1,2))==min(abs(Y-fix(1,2))));
E_1=P(row,col)
row1=find(abs(X-fix(3,1))==min(abs(X-fix(3,1))));
col1=find(abs(Y-fix(3,2))==min(abs(Y-fix(3,2))));
E_2=P(row1,col1)
E_w=min(min(P))

    s1=find(abs(X-fix(2,1))==min(abs(X-fix(2,1))));
    s2=find(abs(Y-fix(2,2))==min(abs(Y-fix(2,2))));
    E_s=P(s1,s2)
    E_0=P(51,51)
%     save('Ecp-8','E_0','E_1','E_2','E_s','E_w')

sum(sum(Z))