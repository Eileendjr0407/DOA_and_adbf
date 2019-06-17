clear all£»clc;close all;
M=80;  % The number of arrays
theta=-90:1:90;
theta0=40;
a0=exp(j*pi*sind(theta0)*(0:M-1)');
a=exp(j*pi*(0:M-1)'*sind(theta));
p_=conj(a0')*a;
p_=10*log10(abs(p_)/max(abs(p_)));
figure(1)
plot(theta,p_);
grid on