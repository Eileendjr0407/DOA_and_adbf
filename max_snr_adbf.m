%% Title: ADBF(����Ӧ�����γ�) ���õ�׼������������׼��
%% Author�� ������
%%  ����Ӧ�����γ��㷨
%% Date: 2019-6-17
clear all;close all;
% clc
M=64;   % The number of signal
Ns=256 ;   % ��������
Nj=2;
lambda_d=1/2;   % ��Ԫ��ࣨ�����ı�ֵ��
thetaj=[-30,20];   % ���ŷ���
theta0=10;      %  ����ָ��
theta=-90:1:90;   %  ������������
SNR=10;
JNR=60;
j=sqrt(-1);
%% Establish the signal model
nj=length(thetaj);
K=length(theta0);
a0=exp(j*2*pi*lambda_d*[0:M-1]'*sind(theta0));
aj=exp(j*2*pi*lambda_d*[0:M-1]'*sind(thetaj));
a=exp(j*2*pi*lambda_d*[0:M-1]'*sind(theta));
signal=10^(JNR/20)*0.707*(randn(nj,Ns)+j*randn(nj,Ns));
signal_1=10^(JNR/20)*0.707*(randn(K,Ns)+j*randn(K,Ns));
noise=0.707*(randn(M,Ns)+j*randn(M,Ns));
X=aj*signal+noise;
signal_zhen=a0*signal_1;
%% 
Rs=1/Ns*signal_zhen*signal_zhen';
Rn=1/Ns*X*X';
[V,D]=eig(Rs,Rn);
[c,b]=max(diag(D));
Wopt=V(:,b);
pattern2=abs(Wopt'*a);
pattern2=20*log10(pattern2/max(pattern2));
pattern0=abs(a0'*a)+0.000000001;
pattern0=20*log10(pattern0/max(pattern0));
figure(1)
plot(theta,pattern0,'--r',theta,pattern2);
grid on
legend('��ͨ�����γ�','max-snr�����γ�');
xlabel('�Ƕȣ��ȣ�');
ylabel('����');
title('�����γ�(max-snr)');


