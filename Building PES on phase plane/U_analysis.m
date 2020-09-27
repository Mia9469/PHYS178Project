clear;close all;

%1:mu0�Ŷ�����/2:D�Ŷ�����/3:C'�Ŷ�����
c=3;

switch c
%% MU0
    case 1
% load data
for i=1:13
    load(['E-' num2str(floor(i*2.5-2.5)) '.mat']);
    U(i).E=E;
    U(i).E_w=E_w;
    U(i).E_s=E_s;
    if i<6
        U(i).E_c=E_c;
        U(i).E_lowc=E_lowc;        
    end
end


% all 2_sta potential
p=zeros(1,13);p_s=zeros(1,13);
for i=1:13
    p(i)=U(i).E_w;
    p_s(i)=U(i).E_s-U(i).E;
end
figure;
plot(0:2.5:30,p,'*-');
xlabel('Mu_0/Hz')
ylabel('U')
title('ѡ����̬����-���Դ̼�Mu_0')
figure;
plot(0:2.5:30,p_s,'*-');
xlabel('Mu_0/Hz')
ylabel('U')
title('���ݸ߶�-���Դ̼�Mu_0')

% 3_sta
p_d=zeros(1,5);p_c=zeros(1,5);
for i=1:5
    p_d(i)=U(i).E_s-U(i).E_lowc;
    p_c(i)=U(i).E_lowc-U(i).E;
end
figure;
plot(0:2.5:10,p_c,'*-');
xlabel('Mu_0/Hz')
ylabel('U')
title('�м���̬��ѡ����̬�Ʋ�-���Դ̼�Mu_0')
figure;
plot(0:2.5:10,p_d,'*-');
xlabel('Mu_0/Hz')
ylabel('U')
title('�м���̬�����ర���Ʋ�-���Դ̼�Mu_0')

%% D
    case 2
for i=1:10
    load(['Ed-' num2str(i) '.mat']);
    U(i).E=E;
    U(i).E_w=E_w;
    U(i).E_s=E_s;
end

p=zeros(1,10);p_s=zeros(1,10);
for i=1:10
    p(i)=U(i).E_w;
    p_s(i)=U(i).E_s-U(i).E;
end
figure;
plot(0.1:0.1:1,p,'*-');
xlabel('D')
ylabel('U')
title('ѡ����̬����-��ɢϵ��D')
figure;
plot(0.1:0.1:1,p_s,'*-');
xlabel('D')
ylabel('U')
title('���ݸ߶�-��ɢϵ��D')

%% C'
    case 3
for i=1:9
    load(['Ecp-' num2str(i-1) '.mat']);
    U(i).E_w=E_w;
    U(i).E_1=E_1;
    U(i).E_2=E_2;
    U(i).E_s=E_s;
    U(i).E_0=E_0;
end

p_1=zeros(1,9);p_2=zeros(1,9);p_s=zeros(1,9);p_0=zeros(1,9);
for i=1:9
    p_1(i)=U(i).E_1;
    p_2(i)=U(i).E_2;
    p_s(i)=U(i).E_s-U(i).E_w;
    p_0(i)=U(i).E_0;
end
figure;
plot(0:6.4:51.2,p_1,'r*-');hold on;
plot(0:6.4:51.2,p_2,'b*-');hold on;
plot(0:6.4:51.2,p_2-p_1,'k*-');
xlabel('Cprime')
ylabel('U')
title('ѡ����̬����-ƫ���Cprime')
legend({'��ȷ����U_r','��������̬����U_w','���ܲ�U_w-U_r'})
figure;
plot(0:6.4:51.2,p_0,'*-');
xlabel('Cprime')
ylabel('U')
title('��ʼ0������-ƫ���Cprime')
% figure;
% plot(0:6.4:51.21,p_s,'*-');
% xlabel('D')
% ylabel('U')
% title('���ݸ߶�-ƫ���Cprime')
end

