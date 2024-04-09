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
    i_peak(i)=2*R*P_avg(i)/0.5;
    % Peak Electrical amplitude
    Ep(i)=i_peak(i)^2*Tb*0.5;
    % Peak energy (Energy per bit is Ep/2)
    sgma(i)=sqrt(N0/2/Tsamp);
    % noise variance (power spectral density related to sampling time)
    % sgma(i)=i_peak(i)/sqrt(2)*sqrt(nsamp/(2*SNR(i)));
    pt=ones(1,nsamp/2)*i_peak(i);
    pt((nsamp/2+1):nsamp)=0;
    % transmitter filter
    rt=pt;
    % receiver filter matched to pt
    ooK1=round(rand(1,sig_length));
    ooK2=zeros(1,sig_length);
    ook=[ooK1;ooK2];
    OOK=reshape(ook,[1 sig_length*2]);
    % OOK random signal generation （此处产生 OOK 信号：OOK）
    Tx_signal=rectpulse(OOK,nsamp/2)*i_peak(i);
    % Pulse shaping function (rectangular pulse)
    Rx_signal=Tx_signal+normrnd(0,sgma(i),1,sig_length*nsamp);
    % received signal (此处编写在 AWGN 信道下接收信号代码 
    % Rx_signal：接收信号=发送信号+噪声)
    MF_out=conv(Rx_signal,rt)*Tsamp;
    % matched filter output
    MF_out_downsamp=MF_out(nsamp/2:nsamp:end);
    % sampling at end of bit period
    MF_out_downsamp=MF_out_downsamp(1:sig_length);
    % truncation
    Eb=Ep(i)/2;
    % thresholding and demodulation (此处设置一个判决门限来解码 1 和 0，
    % 门限为每比特平均能量。)
    error=0;
    for K=1:sig_length
        if (MF_out_downsamp(K)>=Eb)&(ooK1(K)==0)
            error=error+1;
        end
        if (MF_out_downsamp(K)<Eb)&(ooK1(K)==1)
            error=error+1;
        end
    end
    ber(i)=error/sig_length;
    % bit error calculation（此处计算系统仿真数值误码率：ber）
end

figure;
semilogy(EbN0,ber,'b-*','linewidth',2);
hold on
P=qfunc(sqrt(SNR/0.5));
semilogy(EbN0,P,'r-*','linewidth',2);
% analytical performance, Q-function（此处计算系统理论解析误码率，并画图）
grid on
legend('数值解','解析解');
xlabel('Eb/No, dB');
ylabel('BER');
title('OOK-RZ 调制 BER 曲线');
