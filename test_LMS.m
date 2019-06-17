clear all;clc;close all;
g=100;           %仿真统计次数
N=1024;          %输入信号的序列长度
k=128;           %FIR滤波器长度
pp=zeros(g,N-k);   %将每次循环仿真的误差信号结果存于矩阵PP中，以便求取统计平均。
u=1/256;          %步长因子设为1/256
snr=[3 -3];         %存放输入信号信噪比参数
%生成正弦信号序列
t=1:N;
s=sin(0.5*pi*t);  %生成正弦波信号
xn=zeros(1,N);   %存放输入信号
y=zeros(1,N);    %存放输出信号
w=zeros(1,k);    %存放权值信号
e=zeros(1,N);    %存放误差信号
for type=1:4
    for q=1:g
        noise=rand(1,length(s));
        if type==1
           SNR=snr(1);d=s;
        elseif type==2
           SNR=snr(1); d=s;
        elseif type==3
           SNR=snr(2); d=sqrt(10^(-SNR/10))*noise;
        else
           SNR=snr(2); d=sqrt(10^(-SNR/10))*noise;
        end           
        xn=sqrt(10^(-SNR/10))*noise+s;
        y(1:k)=xn(1:k);
        %LMS算法
        for i=(k+1):N
            XN=xn((i-k+1):(i));
            y(i)=w*XN';
            e(i)=d(i)-y(i);
            w=w+u*e(i)'*XN;
        end
        pp(q,:)=(e(k+1:N)).^2;%求每次仿真后误差信号的平方值
    end
    figure(1);
    subplot(311);
    plot(s(300:450));  %截取信号的一段进行绘图
    title('信号S时域波形');
    if type==1
        subplot(312); plot(xn(300:450));
        title('信号s加噪声后的时域波形(snr=3dB)');
    elseif type==3
        subplot(313); plot(xn(300:450));
        title('信号s加噪声后的时域波形(snr=-3dB)');
    end
    %求取各次循环仿真的误差统计均值，完成Monte Carlo仿真
    for b=1:N-k
        bi(b)=sum(pp(:,b))/g;
    end
    %绘制自适应滤波后的输出信号
    figure(2)
    if type==1
        subplot(311);
        plot(y(300:450));title('自适应滤波后的输出时域信号(snr=3dB,期望信号为正弦信号)');
    elseif type==2
        subplot(312);
        plot(y(300:450));title('自适应滤波后的输出时域信号(snr=3dB,期望信号为正弦信号)');
    elseif type==4
        subplot(313);y=xn-y;%由于期望信号为噪声信号，系统相当于干扰抵消系统
        plot(y(300:450));title('自适应滤波后的输出时域信号(snr=-3dB,期望信号为噪声信号)');
    end
    %绘制误差信号图
    figure(3)
    if type==1
        subplot(311);
        plot(bi(1:100));title('误差均方信号(snr=3dB,期望信号为正弦信号)');
    elseif type==2
        subplot(312);
        plot(bi(1:100));title('误差均方信号(snr=3dB,期望信号为正弦信号)');
    elseif type==3
        subplot(313);
        plot(bi(1:100));title('误差均方信号(snr=-3dB,期望信号为噪声信号)');
    end
end
