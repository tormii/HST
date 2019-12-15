function [Tx,t,f,xMean,GroupDelay] = HST(x , fs,  WindowOpt, Parameter, Mode)
%	------------------- Horizontal Synchrosqueezing Transform -------------------- 
% Authors: Xiaotong Tu and Fucai Li
% email:tormiier@gmail.com,tormii@sjtu.edu.cn;
% https://www.researchgate.net/profile/Xiaotong_Tu2
%Input:
%       
%       x:imput signal
%       fs:sampling frequency/采样频率（Hz）
%       WindowOpt:window function/窗函数选择参数
%           WindowOpt.s：(0.01) 窗函数初始尺度
%           WindowOpt.f0：(0) 窗函数初始中心频率
%           WindowOpt.type：(gauss) 窗函数类型
%       Parameter:频率选择参数
%           Parameter.L：(200) 频率划分个数
%           Parameter.fmin：(最小分辨率) 分析最小频率
%           Parameter.fmax：(奈奎斯特频率) 分析最大频率
%       Mode:(1Ord，2Ord)
%Output:
%       Tx:TFR
%       t:time
%       f:frequency/压缩后的频率（Hz）
%       GroupDelay:group delay;
%---------------------------------------------------------------------------------
% When using this code, please cite our papers:
% Xiaotong Tu, Zhoujie He, Yue Hu, Saqlain Abbas and Fucai Li, Horizontal Synchrosqueezing Transform: Algorithm and Applications, IEEE Sensors Journal
% Author: Xiaotong Tu（Apr.,2019）
%---------------------------------------------------------------------------------
%% 预处理信号
    N = length(x);
%% 参数赋值
    s = WindowOpt.s; type = WindowOpt.type;
    L = Parameter.L; fmin = Parameter.fmin; fmax = Parameter.fmax;
    gamma = sqrt(eps); 
%% SST计算
    %STFT计算
    [Wx,t,f,xMean] = stft(x, fs, WindowOpt, Parameter, 'normal');
    %瞬时频率计算/calculate GD
    if strcmp(Mode, '1Ord')
        WindowOpt.type = '1ord(w)_gauss';
        [dWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        GroupDelay = real(dWx./Wx/(1i));
        for ptr = 1:N
            GroupDelay(:,ptr) = GroupDelay(:,ptr) + t(ptr);
        end
        GroupDelay( abs(Wx) < gamma ) = Inf;
    elseif strcmp(Mode, '2Ord')
        WindowOpt.type = '1ord(w)_gauss';
        [dWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = '2ord(w)_gauss';
%         [ddWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w*gauss';
        [wWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'w*1ord(w)_gauss';
        [wdWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        WindowOpt.type = 'ww*gauss';
        [wwWx,~,~,~] = stft(x, fs, WindowOpt, Parameter, 'normal');
        Denominator = wwWx.*Wx-wWx.*wWx;
        Numerator = dWx.*wWx-Wx.*Wx-Wx.*wdWx;
        q = Numerator./Denominator/(1i);%q = Numerator./Denominator/(1i*2*pi)
        tao_temp=dWx./(1i*Wx);%tao_temp=dWx./(1i*2*pi*Wx)
        for ptr = 1:N
            tao(:,ptr) = -tao_temp(:,ptr) - t(ptr);
        end
        omega_temp=wWx./Wx;
%          for pwr = 1:L
%            omega(pwr,:) = omega_temp(pwr,:) + f(pwr);
%         end
      
        GroupDelay = -real(tao+q.*(-omega_temp));
        GroupDelay( abs(Denominator) < gamma ) = Inf;
    else
        error('Unknown SST Mode: %s', Mode);
    end
    %频率差分计算
    dt = 1/fs;
    %窗时域生成
    [gf,~] = windowf(s,type);
    %计算g(0)
    g0 = gf(0);
    g0 = conj(g0);
    if(g0 == 0)
        error('window must be non-zero and continuous at 0 !');
    end
    %同步压缩/Reassignment
    Wx(isinf(GroupDelay)) = 0;
    Tx = zeros(L,N);
    for prt=1:L
        for b=1:N
            m = min(max(1 + round((GroupDelay(prt,b)-0)/dt),1),N);
            Tx(prt, m) = Tx(prt, m) + Wx(prt, b)*dt;
        end
    end
    Tx = Tx / g0;
end
