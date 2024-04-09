clear all;
close all;
clc;

N=64;                   % Number of Subcarriers
N_data_symbol=N/2-1;    % Number of effective subcarriers for information data
CP=N/4;                 % Length of CP
M=16;                   % Order of modulation
BitperSymbol=log2(M);   % Number of bits in a symbol
N_Iteration=10;       % Number of iterations
SNR_dB=0:1:15;          % Sianal to Noise Ratio（SNR) in dB
Nsym=256;               % Number of DCO-OFDM blocks within a frame
Es=1;                   % Symbol energy is normalized to 1
Eb=Es/BitperSymbol;     % Energy per bit
K=3.2;                  % DC bias ratio

for snr=0:1:15
    snr
    Eb_N0=10.^(snr/10);
    N0=Eb./Eb_N0;
    for k=1:N_Iteration
        %% Input bit streams for modulation
            % 步骤一:编写子函数，将产生的二进制比特数调制成M-QAM码
            [bit,QAM] = Modulation_bit(N_data_symbol,BitperSymbol,N,Nsym);
        %% Hermitian Symmetry
            % 步骤二:编写子函数，将QAM信号变换为具有赫密特对称的形式
            [x]= Hermitain_sym(Nsym,bit);
        %% IFFT
            % 步骤三:编写子函数，将具有赫密特对称的复数信号进行IFFT后，产生双极性时域信号
            [x_ifft] = do_IFFT(x,N);
        %% add DC and clipping
            % 步骤四:编写子函数，产生直流偏置,把直流偏置叠加到时域信号，并截取负信号，让剩余负信号变成0
            [X] = DC_Clip(x_ifft,K);
        %% Guard Interval ingertion and the CP addition
            % 步骤五:编写子函数,时域信号增加循环前缀
            [X1] = add_cp(X,CP,N);
        %% Received sign with AWGN
            % 步骤六:编写子函数，考虑AWGN信道，接收信号只受到噪声影响
            [Y] = Channel(X1,N0);
        %% Remove CP
            % 步骤七:编写子函数,去掉接收信号循环前缀
            [Y] = Remove_CP(Y,N,CP);
        %% FFT
            % 步骤八:编写子函数,接收信号经过FFT后变成频域信号
            [Y1] = do_FFT(Y,N);
        %% demodulation
            % 步骤九:编写子函数,把频域信号在每一个载波上进行解调
            [Y3] = demodulation(Y1,N,Nsym);
        %% BER calculatio
            % 步骤十:编写子函数,计算系统误码率
            [BER(snr+1)] = BER_calculation(Y3,QAM,N,Nsym,snr);
    end
    ber(snr+1)=4*(sqrt(M)-1)*qfunc(sqrt(3*log2(M)*Eb_N0/(M-1)))/(sqrt(M)*log2(M))...
    +4*(sqrt(M)-2)*qfunc(3*sqrt(3*log2(M)*Eb_N0/(M-1)))/(sqrt(M)*log2(M));
end

figure;
semilogy(SNR_dB,BER,'r-*'); %此处添加自己定义的BER变量
hold on
grid on
semilogy(SNR_dB,ber,'b--'); %此处绘制解析BER
legend('实验解','解析解')
xlabel('EbN0(dB)')
ylabel('BER')
title('BER for DCO-OFDM')

function [bit,QAM] = Modulation_bit(N_data_symbol,BitperSymbol,N,Nsym)
    % Modulation_bit 将产生的二进制比特数调制成M-QAM码
    bits_num=N_data_symbol*BitperSymbol*Nsym;
    bits=round(rand(bits_num,1));
    bit=qammod(bits,16,'inputtype','bit');
    bit=bit/sqrt(10);
    bit=reshape(bit,N_data_symbol,Nsym);
    QAM=bits;
end

function [x]= Hermitain_sym(Nsym,bit)
    % Hermitain_sym 将QAM信号变换为具有赫密特对称的形式
    x=[zeros(1,Nsym);bit;zeros(1,Nsym);flipud(conj(bit))];
end

function [x_ifft] = do_IFFT(x,N)
    % do_IFFT 将具有赫密特对称的复数信号进行IFFT后，产生双极性时域信号
    x_ifft=ifft(x)*sqrt(N);
end

function [X] = DC_Clip(x_ifft,K)
    % DC_Clip 产生直流偏置,把直流偏置叠加到时域信号，并截取负信号，让剩余负信号变成0
    X=x_ifft+K;
    X(X<0)=0;
end

function [X1] = add_cp(X,CP,N)
    % add_cp 时域信号增加循环前缀
    X1=[X(N-CP+1:N,:);X];
end

function [Y] = Channel(X1,N0)
    % Channel 考虑AWGN信道，接收信号只受到噪声影响
    Y=X1+sqrt(N0)*randn(size(X1));
end

function [Y] = Remove_CP(Y,N,CP)
    % Remove_CP 去掉接收信号循环前缀
    Y=Y(CP+1:end,:);
end

function [Y1] = do_FFT(Y,N)
    % do_FFT 接收信号经过FFT后变成频域信号
    Y1=fft(Y)/sqrt(N);
    Y1=Y1(2:N/2,:);
end

function [Y3] = demodulation(Y1,N,Nsym)
    % demodulation 把频域信号在每一个载波上进行解调
    Y1 = Y1*sqrt(10);
    Y2 = qamdemod(Y1,16,'OutputType','bit');
    Y3 = reshape(Y2,numel(Y2),1);
end

function [BER] = BER_calculation(Y3,QAM,N,Nsym,snr)
    % BER_calculation 计算系统误码率
    BER=sum(abs(Y3-QAM))/numel(Y3);
end
