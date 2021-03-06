%% Title: ADBF(自适应波束形成) 利用的准则是最小均方误差准则
%% Author： 丁杰如
%%  自适应波束形成算法
%% Date: 2019-6-17
clear all;close all;
% clc
M=64;   % The number of signal
Ns=256 ;   % 采样点数
Nj=2;
K=2;    %  辅助天线阵元数目
lambda_d=1/2;   % 阵元间距（波长的比值）
thetaj=[-30,20];   % 干扰方向
theta0=10;      %  波束指向
theta=-90:1:90;   %  波束搜索区间
SNR=10;
JNR=60;
j=sqrt(-1);
%% Establish the signal model
nj=length(thetaj);
Vs=exp(j*2*pi*lambda_d*(0:M-1)'*sind(theta));
Vs0=exp(j*2*pi*lambda_d*(0:M-1)'*sind(theta0));
Vsj=exp(j*2*pi*lambda_d*(0:M-1)'*sind(thetaj));
AJ=10^(JNR/20)*0.707*(randn(nj,Ns)+j*randn(nj,Ns));
noise=0.707*(randn(M,Ns)+j*randn(M,Ns));
Xs=Vsj*AJ+noise;
Xj=Xs(1:K,:);
%% 
D=Vs0'*Xs;
R11=Xj*Xj'/Ns;
rxd=Xj*D'/Ns;
W=inv(R11)*rxd;
pattern1=abs(Vs0'*Vs-W'*Vs(1:K,:))+0.00000001;    %  对消 最优的滤波器-辅助天线的
% pattern1=abs(W'*Vs)+0.00000001;
pattern1=20*log10(pattern1/max(pattern1));
pattern0=abs(Vs0'*Vs)+0.00000001;
pattern0=20*log10(pattern0/max(pattern0));
figure(1)
plot(theta,pattern0,'--r',theta,pattern1)
legend('普通波束形成','自适应波束形成');
xlabel('角度（度）');
ylabel('幅度');
grid on
title('自适应波束形成（MSE准则）')
%% 