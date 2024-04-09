clear all
close all

q=1.6e-19;
% Charge of Electron
Ib=202e-6;
% Background Noise Current + interfernce
N0=2*q*Ib;
% Noise Spectral Density （此处需要产生噪声功率谱密度：N0）
R=1;
% Photodetector responsivity
Rb=1e6;
% Bit rate
Tb=1/Rb;
% bit duration（此处需要计算比特时间：Tb）
sig_length=1e5;
% number of bits
nsamp=10;
% samples per symbol
Tsamp=Tb/nsamp;
% sampling time
EbN0=1:12;
% signal-to-noise ratio in dB.
SNR=10.^(EbN0./10);
% signal-to-noise ratio

% ********** Simulation of probability of errors. ************
for i=1:length(SNR)
    P_avg(i)=sqrt(N0*Rb*SNR(i)/(2*R^2));
    % average transmitted optical power
    i_peak(i)=2*R*P_avg(i);
    % Peak Electrical amplitude
    Ep(i)=i_peak(i)^2*Tb;
    % Peak energy (Energy per bit is Ep/2)
    sgma(i)=sqrt(N0/2/Tsamp);
    % noise variance (power spectral density related to sampling time)
    %sgma(i)=i_peak(i)/sqrt(2)*sqrt(nsamp/(2*SNR(i)));
    pt=ones(1,nsamp)*i_peak(i);
    % transmitter filter
    rt=pt;
    % receiver filter matched to pt
    OOK=round(rand(1,sig_length));
    % OOK random signal generation （此处产生 OOK 信号：OOK）
    Tx_signal = reshape((pt.'*OOK),1,sig_length*nsamp);
    % Pulse shaping function (rectangular pulse)
    Rx_signal=Tx_signal+normrnd(0,sgma(i),1,sig_length*10);
    % received signal (此处编写在 AWGN 信道下接收信号代码 
    % Rx_signal：接收信号=发送信号+噪声)
    MF_out=conv(Rx_signal,rt)*Tsamp;
    % matched filter output
    MF_out_downsamp=MF_out(nsamp:nsamp:end);
    % sampling at end of bit period
    MF_out_downsamp=MF_out_downsamp(1:sig_length);
    % truncation
    Eb=Ep(i)/2;
    % thresholding and demodulation (此处设置一个判决门限来解码 1 和 0，
    % 门限为每比特平均能量。)
    error1=sum(((MF_out_downsamp-Eb).*OOK)<0);
    error0=sum(((MF_out_downsamp-Eb).*(OOK-1))<0);
    ber(i)=(error0+error1)/sig_length;
    P(i)=(erfc(sqrt(Eb/N0)/sqrt(2)))/2;
    % bit error calculation（此处计算系统仿真数值误码率：ber）
end

figure;
semilogy(EbN0,ber,'b-*','linewidth',2);
hold on
semilogy(EbN0,P,'r-*','linewidth',2);
% analytical performance, Q-function（此处计算系统理论解析误码率，并画图）
grid on
legend('数值解','解析解');
xlabel('Eb/No, dB');
ylabel('BER');
title('OOK-NRZ 调制 BER 曲线');
