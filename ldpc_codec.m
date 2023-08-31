clear all
clc

% Размер блока кодирования (максимальный 98)
N = 68;

% Количество итераций при декодировании LDPC
maxnumiter = 12;

% Позиционность PSK модуляции
M = 8;

% Скорость помехоустойчивого кодирования
rate = 1/2;

% Матрица проверки четности LDPC для скорости кодирования 1/2
parity_check_1_2 = [ -1, 94, 73, -1, -1, -1, -1, -1, 55, 83, -1, -1, 7, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1;
                                -1, 27, -1, -1, -1, 22, 79, 9, -1, -1, -1, 12, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1;
                                -1, -1, -1, 24, 22, 81, -1, 33, -1, -1, -1, 0, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1;
                                61, -1, 47, -1, -1, -1, -1, -1, 65, 25, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1, -1;
                                -1, -1, 39, -1, -1, -1, 84, -1, -1, 41, 72, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1, -1;
                                -1, -1, -1, -1, 46, 40, -1, 82, -1, -1, -1, 79, 0, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1, -1;
                                -1, -1, 95, 53, -1, -1, -1, -1, -1, 14, 18, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1, -1;
                                -1, 11, 73, -1, -1, -1, 2, -1, -1, 47, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1, -1;
                                12, -1, -1, -1, 83, 24, -1, 43, -1, -1, -1, 51, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1, -1;
                                -1, -1, -1, -1, -1, 94, -1, 59, -1, -1, 70, 72, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, -1;
                                -1, -1, 7, 65, -1, -1, -1, -1, 39, 49, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0;
                                43, -1, -1, -1, -1, 66, -1, 41, -1, -1, -1, 26, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0];

% Матрица проверки четности LDPC для скорости кодирования 3/4
parity_check_3_4 = [6, 38, 3, 93, -1, -1, -1, 30, 70, -1, 86, -1, 37, 38, 4, 11, -1, 46, 48, 0, -1, -1, -1, -1;
                                62, 94, 19, 84, -1, 92, 78, -1, 15, -1, -1, 92, -1, 45, 24, 32, 30, -1, -1, 0, 0, -1, -1, -1;
                                71, -1, 55, -1, 12, 66, 45, 79, -1, 78, -1, -1, 10, -1, 22, 55, 70, 82, -1, -1, 0, 0, -1, -1;
                                38, 61, -1, 66, 9, 73, 47, 64, -1, 39, 61, 43, -1, -1, -1, -1, 95, 32, 0, -1, -1, 0, 0, -1;
                                -1, -1, -1, -1, 32, 52, 55, 80, 95, 22, 6, 51, 24, 90, 44, 20, -1, -1, -1, -1, -1, -1, 0, 0;
                                -1, 63, 31, 88, 20, -1, -1, -1, 6, 40, 56, 16, 71, 53, -1, -1, 27, 26, 48, -1, -1, -1, -1, 0];

% Приведение матрицы кодирования к размеру блока
if rate == 1/2
   parity_check = round(parity_check_1_2./95.*N);
elseif rate == 3/4
    parity_check = round(parity_check_3_4./95.*N);
else
    parity_check = round(parity_check_1_2./95.*N);
end

ber = comm.ErrorRate;

% Формирование проверочной матрицы
v = ones(1, N);
block = diag(v, 0);
block_zero = zeros(N,N);
parity_check_matrix = [];
for row = 1 : size(parity_check, 1)
    prom = [];
    for col = 1 : size(parity_check, 2)
        ind = parity_check(row, col);
        if ind == -1
            prom = [prom block_zero];
        else
            Y = circshift(block,parity_check(row, col));
            prom = [prom Y];
        end
    end
    parity_check_matrix = [parity_check_matrix; prom];
end

% Приведение матрицы к логическому виду
rez = logical(sparse(parity_check_matrix));

% Формируем кодер и декодер LDPC(1632, 816)
cfgLDPCEnc = ldpcEncoderConfig(rez);
cfgLDPCDec = ldpcDecoderConfig(rez);

% Параметры БЧХ кода
Nbch=204;
Kbch=188;
% Формируем кодер и декодер БЧХ(204, 188)
bchEncoder = comm.BCHEncoder(Nbch, Kbch);
bchDecoder = comm.BCHDecoder(Nbch, Kbch);

rezSNR = [];
ind = 1;
for SNR = 0 : 0.2 : 5
    SNR
    errors = 0;
    info_bits = 0;
    for num_code_word = 1 : 1000
        % Формирование информационных данных
        if rate == 1/2
            dataOut = randi([0 1], Kbch*4,1, 'int8');
        elseif rate == 3/4
            dataOut = randi([0 1], Kbch*6,1, 'int8');
        else
            dataOut = randi([0 1], Kbch*4,1, 'int8');
        end

        % 1-ый каскад. Кодирование кодом БЧХ
        dataEncodedBCH = bchEncoder(dataOut);
        
        % 2-ой каскад. Кодирование LDPC
        encodedData = ldpcEncode(dataEncodedBCH,cfgLDPCEnc);
        
        % Модуляция сигнала
        modSignal = pskmod(encodedData,M,InputType='bit');
        
        % Добавляем к сигналу шум
        [rxsig, noisevar] = awgn(modSignal, SNR, 'measured');
        
        % Демодуляция сигнала
        demodSignal = pskdemod(rxsig, M, OutputType='approxllr', NoiseVariance=noisevar);
        
        % Снятие кодера LDPC
        dataDecodeLDPC = ldpcDecode(demodSignal,cfgLDPCDec, maxnumiter);
        
        % Снятие кодера BCH
        dataDecodeBCH = bchDecoder(dataDecodeLDPC);
         
        % Получить количество ошибок
        errors = errors + sum(abs(dataOut - dataDecodeBCH));
        % Получить количество информационных бит
        info_bits = info_bits + length(dataOut);
    end

    % Формирование массива SNR/ошибка
    if errors == 0
        rezSNR(ind, :) = [SNR log10(1e-10)];
    else
        rezSNR(ind, :) = [SNR log10(errors/info_bits)];
    end

    ind = ind + 1;
end

figure(100)
plot(rezSNR(:,1), rezSNR(:,2))
title('PSK8 3/4 (RS(204,188), LDPC(1632,1224))')
xlabel('SNR')
ylabel('BER 10^x')
grid on
