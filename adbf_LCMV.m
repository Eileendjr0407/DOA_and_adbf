%% Title: ADBF(����Ӧ�����γ�) ���õ�׼����LCMV(linear constrained Minimum variable)
%% Author�� ������
%%  �����γ��㷨
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
a0=exp(j*2*pi*lambda_d*[0:M-1]'*sind(theta0));
aj=exp(j*2*pi*lambda_d*[0:M-1]'*sind(thetaj));
a=exp(j*2*pi*lambda_d*[0:M-1]'*sind(theta));
signal=10^(JNR/20)*0.707*(randn(nj,Ns)+j*randn(nj,Ns));
noise=0.707*(randn(M,Ns)+j*randn(M,Ns));
X=aj*signal+noise;
%% 
Rx=X*X'*1/Ns;
Wopt=inv(Rx)*a0*inv(a0'*inv(Rx)*a0);
pattern0=abs(a0'*a+0.00000001);
pattern1=abs(Wopt'*a);
pattern0=20*log10(pattern0/max(pattern0));
pattern1=20*log10(pattern1/max(pattern1));
figure(1)
plot(theta,pattern0,'--r',theta,pattern1);
grid on
legend('��ͨ�����γ�','LCMV�����γ�');
xlabel('�Ƕȣ��ȣ�');
ylabel('����');
title('�����γ�(LCMV)');