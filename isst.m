function b = isst(Tx, xMean) 
%	------------------- 短时傅里叶逆变换 -------------------- 
%Input:
%       Wx:短时傅里叶变换系数
%       xMean:正变换时去掉的平均值
%       fs:采样频率（Hz）
%       WindowOpt:窗函数选择参数
%           WindowOpt.s：(0.01) 窗函数初始尺度
%           WindowOpt.f0：(0) 窗函数初始中心频率
%           WindowOpt.type：(gauss) 窗函数类型
%       Parameter:频率选择参数
%           Parameter.L：(200) 频率划分个数
%           Parameter.fmin：(最小分辨率) 分析最小频率
%           Parameter.fmax：(奈奎斯特频率) 分析最大频率
%Output:
%       b:重构信号
%---------------------------------------------------------------------------------
%    Synchrosqueezing Toolbox
%    Authors: 何周杰（2019/2/8）
%---------------------------------------------------------------------------------
%% 逆变换
    %重构
    b = sum(real(Tx),1);
    % 添加平均值
    b = b+xMean;
end