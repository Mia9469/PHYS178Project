close all;clear
load("trial1.mat");
load("trial2.mat");
m=figure;
% for i =1:3
%     if (i==1 || i==3)
%         load('hist_d30_r.mat')
%     else
%         load('hist_s30_r.mat')
%     end
%     z=ones(size(hist))*(i-1);
%     mesh(edges_r(1,2:end),edges_r(2,2:end),z,hist_r);hold on;
%     zlim([0,2]);
%     xlabel("S2");ylabel("S1");
% end

for i =1:3
    if (i==1 || i==3)
        load('3_s.mat')
    else
        load('2_s.mat')
        P=(P+5.7674)/3-4.634;
    end
    z=ones(size(P))*(i-1)*2;
    mesh(X,Y,z,P');hold on;
    
    xlabel("S2");ylabel("S1");zlabel("t/s");
end

h=0.5;
Z=linspace(-h,4,length(S1));

load("trial1.mat");
plot3(S1,S2,Z,'r-','LineWidth',2);hold on;
plot3(S1(1),S2(1),-h,'y*','LineWidth',3);hold on;

load("trial2.mat");
plot3(S1,S2,Z,'g-','LineWidth',2);hold on;
plot3(S1(1),S2(1),-h,'y*','LineWidth',3);hold on;
set(gca,'ZDir','reverse')





