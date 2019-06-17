%% ESPRIT ALGORITHM USED FOR DOA ESTIMATION
%% AUTHOR :丁杰如(特征值分解)
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
% SteerringVect_1=exp(j*2*pi*(0:2:M-1)'*sind(doa)/lambda);   % M*K
% SteerringVect_2=exp(j*2*pi*(1:2:M-1)'*sind(doa)/lambda);   % M*K
% SteerringVect=[SteerringVect_1;SteerringVect_2];
SteerringVect=exp(j*2*pi*(0:M-1)'*sind(doa)*delta_d/lambda);
signal=exp(j*w*(1:T));
Y=SteerringVect*signal;
Y=Y+awgn(Y,SNR);    % The observed signal;
%% The ESPRIT ALGORITHM
X_1=Y(1:sub_M,:);
X_2=Y(2:M,:);
X=[X_1;X_2];
R=X*(X')/T;
[U,S,V]=svd(R);
Us=U(:,1:K);
Us_1=Us(1:sub_M,:);
Us_2=Us(sub_M+1:sub_M*2,:)
M=pinv(Us_1)*Us_2;
[Vm,Dm]=eig(M);
disp(Dm);
estimated_source_doa=-asin(angle(diag(Dm'))/pi)*180/pi;
disp(estimated_source_doa);