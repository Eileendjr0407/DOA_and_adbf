% TLS_ESPRIT ALOGRITHM
% DOA ESTIMATION BY TLS_ESPRIT ALOGRITHM
clear variables;
close all;  clc;
%信噪比(dB)
snr=-10:2:10;
lsnr = length(snr);
RMSE=zeros(1,lsnr);
test_num = 300;

source_number=1; % 信源数
sensor_number=8; % 原阵元数
m=4;             % 子阵元数
N_x=1024;        % 信号长度
snapshot_number=N_x; % 快拍数
w=pi/4;              % 信号频率
l=2*pi*3e8/w;        % 信号波长  
d=0.5*l;             % 阵元间距
array_distance=d;    % 两组阵元轴线方向的间距

source_a= 50;  % 两个信号的入射角度
source_doa= source_a *pi/180;  % 两个信号的入射角度
% 阵列流型
A=[exp(-1i*(0:sensor_number-1)*d*2*pi*sin(source_doa)/l)].';

for isnr = 1:lsnr
    RMSE_i = zeros(1, test_num);
    for iter = 1:test_num
        % 仿真信号
        s=sqrt(10.^(snr(isnr)/10))*exp(1i*w*[0:N_x-1]); 
        %x=awgn(s,snr);
        % 加了高斯白噪声后的阵列接收信号
        x=A*s+(1/sqrt(2))*(randn(sensor_number,N_x)+1i*randn(sensor_number,N_x));
        x1=x(1:m,:);     % 子阵1接受的数据矢量
        x2=x(2:(m+1),:); % 子阵2接受的数据矢量

        % 对两个子阵的模型进行合并
        X = [x1; x2];
        R=X*X'/snapshot_number;
        % 对R进行奇异值分解
        [U,S,V]=svd(R);
        R=R-S(2*m,2*m)*eye(2*m);
        [U,S,V]=svd(R);
        Us=U(:,1:source_number);
        % disp(Us);
        Us1=Us(1:m,:);
        Us2=Us((m+1):2*m,:);
        %Us12=[Us1 Us2];
        % 形成矩阵Us12
        Us12=[Us1,Us2];
        % 对“Us12'*Us12”进行特征分解，得到矩阵E
        [E,Sa,Va]=svd(Us12'*Us12);
        % disp('E');
        % disp(E);
        % disp(Sa);
        % 将 E_4x4 分解为四个小矩阵
        E11=E(1,1);
        E12=E(1,2);
        E21=E(2,1);
        E22=E(2,2);
        % 按照公式得到旋转不变矩阵M
        M=-(E12*(inv(E22)));
        % disp('M');
        % disp(M);
        % 对得到的旋转不变矩阵进行特征分解
        [Vm,Dm]=eig(M);
        % disp(Dm);
        Dm=(diag(Dm)).';
        doa=-asin(angle(Dm)/pi)*180/pi;
        disp(doa);
        RMSE_i(iter) = (doa-source_a)^2;
    end
    RMSE(isnr) = sqrt(sum(RMSE_i)/test_num);
end

%% plot 
figure;
plot(snr,RMSE,'m-s');
xlabel('SNR (dB)');
ylabel('RMSE of(degree)');
grid on;  
