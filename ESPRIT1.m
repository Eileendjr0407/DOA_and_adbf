%% ESPRIT ALGORITHM USED FOR DOA ESTIMATION
%% AUTHOR :丁杰如(广义特征值分解)
%% data: 2016-6-16
clear all;
clc; close all;
%% establish the model
M=8;    % The number of arrays
sub_M=M-1;
K=3;    % The number of signal source
c=3e8;    %  the speed of light
Fc=3e9;
j=sqrt(-1);
lambda=c/Fc;   % wavelength
delta_d=lambda/2;   % the distance between two arrays
doa=[10,30,50];     % the angle 
w=[1,2,3]';     % the independent signal source
T=1040;        %  the number of snapshot
SNR=10;      % signal-Noise ratio
SteerringVect=exp(j*2*pi*(0:M-1)'*sind(doa)*delta_d/lambda);
signal=exp(j*w*(1:T));
Y=SteerringVect*signal;
Y=Y+awgn(Y,SNR);    % The observed signal;
X_1=Y(1:sub_M,:);
X_2=Y(2:M,:);
X=[X_1;X_2];
%% 
Rxx=X_1*X_1'/T;
Rxy=X_1*X_2'/T;
% [Vm,Dm]=eig(Rxx);
% data=diag(Dm);
data=eig(Rxx);
Sigma=mean(data(sub_M-K,:));
Cxx=Rxx-Sigma*eye(sub_M);
e=eig(Cxx,Rxy);   % 求广义特征值
% [v,d]=eig(Cxx/Rxx);
% data1=diag(d);
[Y L]=sort(abs(abs(e)-1));
D=e(L);
D_big=D(1:K);
estimated_source_doa=asin(angle((D_big))/pi)*180/pi
