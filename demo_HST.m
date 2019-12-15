  
%%Comparisions/仿真信号TFR不同方法对比
clc; clear all; close all;    
    N = 2000;
    fs = 200;
    t = (0:N-1)/fs;
    f = (0:N/2)*fs/N;
  
%% generate signal/信号生成
%Signal 1/模式1
    a1 = 0;
    b1 = 8;
    c1 = -0.05;
    d1 = 0;
    e1 = 0;
    A_f1 = exp(0.005*f);
    GD_t1 = b1+c1*f+d1*f.^2+e1*f.^3;
    Phi_f1 = a1+b1*f+c1*f.^2/2+d1*f.^3/3+e1*f.^4/4;    
    X1 = A_f1.*exp(-1i*2*pi*Phi_f1);
    X1(end) = -A_f1(end);
    Y1 = [X1  conj(fliplr(X1(2:end-1)))];    
    y1 = ifft(Y1);
%Signal 2/模式2
    a2 = 0;
    b2 = 3;
    c2 = 0.01;
    d2 = 0.003;
    e2 = -3*10^-5;
    A_f2 = exp(0.008*f);
    GD_t2 = b2+c2*f+d2*f.^2+e2*f.^3;
    Phi_f2 = a2+b2*f+c2*f.^2/2+d2*f.^3/3+e2*f.^4/4;    
    X2 = A_f2.*exp(-1i*2*pi*Phi_f2);
    X2(end) = -A_f2(end);%
    Y2 = [X2  conj(fliplr(X2(2:end-1)))];    
    y2 = ifft(Y2);
%Signal 3/模式3
    a3 = 0;
    b3 = 31/10;
    c3 = 7/50;
    d3 = -1/1000;
    e3 =0;
    A_f3 = exp(0.01*f);
    GD_t3 = b3+c3*f+d3*f.^2+e3*f.^3;
    Phi_f3 = a3+b3*f+c3*f.^2/2+d3*f.^3/3+e3*f.^4/4;    
    X3 = A_f3.*exp(-1i*2*pi*Phi_f3);
    X3(end) = -A_f3(end);%
    Y3 = [X3  conj(fliplr(X3(2:end-1)))];    
    y3 = ifft(Y3);
% add noise/加噪声
         y = y2;
%     y = y2;
%      y = awgn(y,noise,'measured');         %为信号加载白噪声
   figure;
   fsize=12;
   set(gcf,'Position',[50 50 560 840]);
   subplot(2,1,1)
    plot(t,y);
    xlabel({'Time (s)','(a)'},'FontSize',10);
    ylabel('Amplitude','FontSize',10);
    set(gca,'fontsize',fsize,'linewidth',1)
    subplot(2,1,2)
    plot(GD_t2,f);
    axis([0 10 0 100])
     xlabel({'Time (s)','(b)'},'FontSize',10);
     ylabel('Frequency(Hz)','FontSize',10);
     set(gca,'fontsize',fsize,'linewidth',1)
%% imput parameters/输入参数
    noise = 0;             %信号高斯噪声分贝，值越小噪声越大
    fs = 200;              %sampling frequency/信号采样率
    N = 2000;               %信号点数
    %SST计算阶数选择
    Mode = '2Ord';        %(1Ord，2Ord)
    %window function/窗选择参数
    WindowOpt = struct('type','gauss','s',0.12);%0.05
    %parameters/频率选择参数
    Parameter = struct('L',N/2+1,'fmin',0,'fmax',fs/2);
%         y = awgn(y,noise,'measured');         %为信号加载白噪声
    %%  STFT 
    [Wx,t1,f1,~] = stft(y, fs, WindowOpt, Parameter, 'normal');
     r = renyi(abs(Wx),t,f',3)
    figure;
    set(gcf,'Position',[50 50 560*2 840]);
    subplot(2,2,1)
    imagesc(t1,f1,abs(Wx));axis xy;
    text(0,90,['Rényi:',num2str(r,3)],'Color','r','FontWeight','bold','FontSize',20);
    title('STFT','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    rectangle('Position',[7.8,58,0.5,20],'EdgeColor','r');

  
    %% TSST
    Mode = '1Ord';
    [hst,t,f,xMean,GD1] = HST(y , fs,  WindowOpt, Parameter, Mode);
    r1 = renyi(abs(hst),t,f',3)
    subplot(2,2,2)
    imagesc(t,f,abs(hst));axis xy;
    text(0,90,['Rényi:',num2str(r1,3)],'Color','r','FontWeight','bold','FontSize',20);
    title('TSST','FontSize',14);
    axis tight; xlabel('Time (t)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    rectangle('Position',[7.8,58,0.5,20],'EdgeColor','r');
   
% %% 参数选择 
%% HST
 Mode = '2Ord';
    [hst,t,f,xMean,GD2] = HST(y , fs,  WindowOpt, Parameter, Mode);
    r2 = renyi(abs(hst),t,f',3)
    subplot(2,2,3)
    imagesc(t,f,abs(hst));axis xy;
    text(0,90,['Rényi:',num2str(r2,3)],'Color','r','FontWeight','bold','FontSize',20);
    title('HST','FontSize',14);
    axis tight; xlabel('Time (t)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    rectangle('Position',[7.8,58,0.5,20],'EdgeColor','r');
    %% SST
    [sst,t,f,xMean,~] = sst(y , fs,  WindowOpt, Parameter, Mode);
     r = renyi(abs(sst),t,f',3)
    subplot(2,2,4)
    imagesc(t,f,abs(sst));axis xy;
      text(0,90,['Rényi:',num2str(r,3)],'Color','r','FontWeight','bold','FontSize',20);
    title('SST','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Frequency(Hz)','FontSize',10);
    rectangle('Position',[7.8,58,0.5,20],'EdgeColor','r');

%%  Reconstruct mode by HST /重构 HST
   singleside=5;
   direction=3;
    Real_GD=GD_t2;
    Ex_HST = ExtractTSST_new( hst,Real_GD,fs,singleside,direction);
%     figure;
%    imagesc(t,f,abs(Ex));axis xy;
%    title('HST-Filtered','FontSize',14);
%    axis tight; xlabel('Time (t)','FontSize',10);
%    ylabel('Frequency(Hz)','FontSize',10);
   
    b = itsst(Ex_HST,fs, xMean);
    figure
     set(gcf,'Position',[50 50 560 840]);
     subplot(2,1,1);
    plot(t,y2,t,b);
    SNRoutput_hst = SNR(real(y2),b)
    
    title('Reconstruction-HST','FontSize',14);
    axis tight; xlabel('Time (s)','FontSize',10);
    ylabel('Amplitude','FontSize',10);
    legend('Origin','New');

   
%% Reconstruct mode by SST重构
    Real_GD=GD_t2;
    Ex_SST = ExtractTSST_new( sst, Real_GD, fs,singleside,direction);
%     figure;
%    imagesc(t,f,abs(Ex));axis xy;
%    title('SST-Filtered','FontSize',14);
%    axis tight; xlabel('Time (t)','FontSize',10);
%    ylabel('Frequency(Hz)','FontSize',10);
   
   b = isst(Ex_SST, xMean);
   SNRoutput_sst = SNR(real(y2),b)
    subplot(2,1,2);
    plot(t,y2,t,b);
    title('Reconstruction-SST','FontSize',14);
    axis tight; xlabel('Time (t)','FontSize',10);
    ylabel('Amplitude','FontSize',10);
    legend('Origin','New');
   

   
    