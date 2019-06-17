%% Title: ADBF(自适应波束形成) 利用的准则是LCMV(linear constrained Minimum variable)
%% Author： 丁杰如
%%  波束形成算法
%% Date: 2019-6-17
clear all;close all;
% clc
M=64;   % The number of signal
Ns=256 ;   % 采样点数
Nj=2;
lambda_d=1/2;   % 阵元间距（波长的比值）
thetaj=[-30,20];   % 干扰方向
theta0=10;      %  波束指向
theta=-90:1:90;   %  波束搜索区间
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
legend('普通波束形成','LCMV波束形成');
xlabel('角度（度）');
ylabel('幅度');
title('波束形成(LCMV)');