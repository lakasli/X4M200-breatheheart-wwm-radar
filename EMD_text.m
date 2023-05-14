clear all; clc; close all;
%O_Signal=load('30s.txt');   % 读入信号 

%雷达模拟参数设置
A = 1;                  %发射信号的振幅
T = 4;              %信号时宽
B = 80;               %信号带宽
F0 = 2*pi;                 %中频频率，即载频频率
N = 2000;
Fs = 1000;       %采样频率
st = LFM_signal(A,T,B,F0);%线性调频信号作为待分析非平稳信号
noise = -5;%高斯白噪声
O_Signal = awgn(st,noise);  %调用 awgn 函数生成一个加入高斯白噪声的信号，并将结果存储在 O_Signal 变量中

figure(2);
subplot(2,1,1);
t=linspace(0,T/2,N);
plot(t,O_Signal);

subplot(2,1,2)
freq=linspace(-Fs/2,Fs/2,N);
plot(freq,fftshift(abs(fft(O_Signal))));

%设定参数
IMF = [];
R_Signal = [];%创建IMF和rt可变列表
i = 0;
sigma1 = 0.1;
sigma2 = 0.1;
alphar1 = 0.1;

%EMD分解
while(1)
    N = length(O_Signal);% 数据长度
    n = 1:N;       % 设置样点序列
    detrend_O_Signal = detrend(O_Signal);
    [O_Signal_MaxLocal,O_Signal_Maxval] = v_findpeaks(detrend_O_Signal,'q',0); % 求极大值位置和幅值
    up_O_Signal = spline(O_Signal_MaxLocal,O_Signal_Maxval,t);         % 内插,获取上包络曲线
    [O_Signal_MinLocal,O_Signal_Minval] = v_findpeaks(detrend_O_Signal,'v',0);% 求极小值位置和幅值
    down_O_Signal = spline(O_Signal_MinLocal,O_Signal_Minval,n);       % 内插,获取下包络曲线
%     K1_length = length(O_Signal_MaxLocal);%获取上包络的极大值个数
%     K2_length = length(O_Signal_MinLocal);%获取上包络的极小值个数
    M_value = (up_O_Signal + down_O_Signal)/2;%m(t)=上包络+下包络的均值
    C_Signal = O_Signal - M_value;%原始信号s（t）-m（t）= c（t）
    figure;
    hold on;
    plot(t,up_O_Signal,'m');
    plot(t,down_O_Signal,'g');
    plot(t,M_value,'r');
    plot(t,O_Signal,'b');
    plot(t,C_Signal,'k');

%     detrend_C_Signal = detrend(C_Signal);
%     [C_Signal_MaxLocal,C_Signal_Maxval] = v_findpeaks(detrend_C_Signal,'q',0); % 求极大值位置和幅值
%     C_Max_value = sum(C_Signal_Maxval) / length(C_Signal_MaxLocal) ;    % 内插,获取上包络曲线均值
%     [C_Signal_MinLocal,C_Signal_Minval] = v_findpeaks(detrend_C_Signal,'v',0);% 求极小值位置和幅值
%     C_Min_value = sum(C_Signal_Minval) / length(C_Signal_MinLocal) ;       % 内插,获取下包络曲线均值


    C_Signal_y = diff(C_Signal);            % 计算一阶差分
    sgn = sign(C_Signal_y);          % 计算差分符号的变化情况
    zero_crossings = find(diff(sgn) ~= 0);  % 找到相邻元素符号不同的位置
    zero_length = length(zero_crossings);  % 统计零点个数
    l = length(findpeaks(C_Signal));
    Diff_length = l - zero_length;
    Ma = abs((up_O_Signal - down_O_Signal))/2;
    p = abs(M_value/Ma);

    idx = find(abs(p) <= sigma1);  % 幅值在 sigma1 以下的位置信息
    durations = diff(idx);  % 计算相邻位置之间的距离
    a_duration = sum(durations);  % 幅值在 sigma1 以下的时间段长度之和
    total_duration = length(C_Signal) - 1;  % 总时间段长度
    ratio = a_duration / total_duration;  % 幅值在 sigma1 以下的时间段占总时间段比值
        
    if ratio>alphar1
        break
    end

%if  (C_Max_value == 0) || ( C_Min_value == 0)  %判断 c(t) 是否满足 IMF 分量的定义
        if (Diff_length==1) || (Diff_length==0)
            i = i+1
            IMF(i) = C_Signal;
            R_Signal(i) = O_Signal - IMF(i);
            O_Signal = R_Signal(i);
            plot(t,O_Signal);
            plot(freq,fftshift(abs(fft(O_Signal))));
        else
            O_Signal = C_Signal;
        end
%     else
%        O_Signal = C_Signal;
%     end
end

figure3;
plot(t,detrend_O_Signal,'k'); 
hold on; 
grid;
plot(t,up_O_Signal,'r');
plot(t,down_O_Signal,'r');
xlabel('样点'); ylabel('幅值');
title('用求取极大极小值方法获取包络曲线图')
set(gcf,'color','w');