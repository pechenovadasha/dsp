function [dX] = FIR_filter_paper(X_t, A, N_used, TH, max_EVM, Nfft, scen)
Nsym=size(X_t,1); % Number of OFDM symbols in the initial time-domain signal

Ndac=scen.Ndac; % Number of digital-to-analog converters
Nant=scen.Ntx; % Number of antennas (> Ndac, HBF architecture)

N_samples=16; % Number of points at the interval 
              % (1 interval - 1 peak per iteration)

coef1=0; % Use only peaks in inversion
coef2=Nant/Ndac; % Use peaks and zeros in inversion
                 % Default ratio for no power loss for coef2: Nant/Ndac 
                 % (Since a peak is a delta function we save only Ndac/Nant
                 % power when cut DAC domain Nant/Ndac times)

Nzero=Nfft-N_used; % padding for IFFT calculation

S_t_ant = permute(X_t, [1 3 2]); % Initial signal to reduce PAPR in
                                 % Size: [Nsym, Nant, Nfft]


% signal power and root mean square
mean_power=mean( abs(S_t_ant(:,:,:)).^2,'all');

% Tunable threasholds (step-by-step reduction can be implemented)
% with increasing the TH, the number of selected peaks is reduced,
% resulting in less PAPR reduction, but less Error vector magnitude (EVM)
%TH = [2.2 1.3];  % SET the THREASHOLD VALUE !!!
TH = [2.5016 2.5119 2.5119];
TH_abs=TH*sqrt( mean_power );
N_iter = length(TH);
% N_iter = 3;

% % % % % % % % % % % % % % Create filter %%%%%%%%%%%%%%%%%%%%
    
filter_order = 32;  % Order of FIR filter
cutoff_frequency = 0.09;  % Example cutoff frequency (0.2 * Nyquist)
fir_coeff = fir1(filter_order, cutoff_frequency);
fir_coeff = fir_coeff/max(abs(fir_coeff));

new_length = 2048;

% Создание нового фильтра с нулями
fir_coeff_padded = zeros(1, new_length);

% Вычисление позиции для размещения исходного фильтра в центре
start_idx = floor((new_length - length(fir_coeff)) / 2) + 1;
end_idx = start_idx + length(fir_coeff) - 1;

% Вставка исходных коэффициентов в центр нового фильтра
fir_coeff_padded(start_idx:end_idx) = fir_coeff;

% Выполнение циклического сдвига
fir_coeff_sinc = circshift(fir_coeff_padded, -new_length / 2);


% % generate SINC signal
% SINC_f=circshift([ones(1,N_used) zeros(1,Nzero)],[0 -(N_used)/2]);
% SINC_t=ifft(SINC_f)*sqrt(Nfft);
% SINC_t=SINC_t/SINC_t(1); % Amplitude normalization to match peak amplitude
% 
% % Generate set of sincs corresponding to different time positions
% % to match the peaks positions
% SINC_mtx = zeros(Nfft, Nfft);
% for j=1:Nfft
%     SINC_mtx(j,:) = circshift(SINC_t,[0 j-1]);
% end

FIR_mtx = zeros(Nfft, Nfft);
for j=1:Nfft
    FIR_mtx(j,:) = circshift(fir_coeff_sinc,[0 j-1]);
end

S_t_canc=zeros(Nsym,Nant,Nfft);
S_t_ant_new=zeros(Nsym,Nant,Nfft);
S_t_ant_canc_peak=zeros(Nsym,Nant,Nfft);

% Number of iterations for peaks positions search per interval
M_iter = 2;

%%%%%%%%%%%%%%New way to find intervals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find peaks to reduce
for i1=1:Nsym
    % fprintf("i1 =%d\n", i1);
    for i2=1:Nant
        % fprintf("i2 =%d\n", i2);
        % Extract signal at particular symbol at particular antenna
        S_t = squeeze(S_t_ant(i1,i2,:)).';
        % Извлекаем интервалы для текущего символа и антенны
        min_inds = find_intervals(S_t);

        % Пошаговое подавление пиков
        for j = 1:N_iter
            S_t_canc_tmp = zeros(1, Nfft);

            % Проходим по всем интервалам
            for k = 1:length(min_inds) - 1
                interval_start = min_inds(k);
                interval_end = min_inds(k + 1) - 1; % Интервалы закрыты

                % Извлекаем сигнал на интервале
                signal = S_t(interval_start:interval_end);

                % Находим максимум и его индекс
                [Max_value, Indx] = max(abs(signal));

                % Проверяем порог
                if Max_value > TH_abs(j)
                    % Амплитуда подавления
                    sinc_Ampl = signal(Indx) * (1 - TH_abs(j) / Max_value);

                    % Глобальная позиция пика
                    sinc_shift = interval_start + Indx - 1;

                    % Суммируем синусоиды для подавления
                    S_t_canc_tmp = S_t_canc_tmp + sinc_Ampl * circshift(fir_coeff_sinc, [0, sinc_shift - 1]);

                    % Запоминаем подавляемый пик
                    S_t_ant_canc_peak(i1, i2, sinc_shift) = ...
                        S_t_ant_canc_peak(i1, i2, sinc_shift) + sinc_Ampl;
                end
            end

            % Обновляем сигнал после подавления
            S_t = S_t - S_t_canc_tmp;
            S_t_canc(i1, i2, :) = squeeze(S_t_canc(i1, i2, :)) + S_t_canc_tmp.';
        end

        % Сохраняем модифицированный сигнал
        S_t_ant_new(i1, i2, :) = S_t;

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % Now we need to check what happens in DAC domain
S_t_dac_canc_sig=zeros(Nsym,Ndac,Nfft);
S_t_dac_canc_peak = zeros(Nsym,Ndac,Nfft);
S_t_ant_canc_peak_new=zeros(Nsym,Nant,Nfft);

for i1=1:Nsym
    for i3=1:Nfft
        % Extract signal at all antennas
        sig = squeeze(S_t_ant_canc_peak(i1,:,i3));

        % Check for found peaks
        indx = find(abs(sig));

        if sum(indx)>0
            A_new = A(indx,:);
            sig_new = sig(indx);
            s1 = pinv(A_new) * sig_new.'; % peaks only transform
            s2 = pinv(A) * sig.'; % peaks + zeros transform
            S_t_dac_canc_peak(i1,:,i3)=coef1*s1+coef2*s2;
        else
            S_t_dac_canc_peak(i1,:,i3) = zeros(1,Ndac);
        end

    end
    % Go back to antenna domain 
    S_t_ant_canc_peak_new(i1,:,:) = A * squeeze(S_t_dac_canc_peak(i1,:,:));

    % Convolve peaks with sinc functions in DAC domain to satisfy spectrum
    % mask
    for i2=1:Ndac
        s_t = squeeze(S_t_dac_canc_peak(i1,i2,:));

        S_t_dac_canc_sig(i1,i2,:) = s_t.' * FIR_mtx;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for EVM and reduce CFR noise proportionally
EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2,'all') / sum(abs(S_t_ant).^2,'all'));
S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx/max_EVM,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the final signal at antennas to analyze PAPR
S_t_ant_new2 = zeros(Nsym,Nant,Nfft);
for i1=1:Nsym
    S_t_ant_new2(i1,:,:) = A * squeeze(S_t_dac_canc_sig(i1,:,:));
end

dS = S_t_ant_new2;
dX = permute(dS, [1,3,2]);

 % test  = S_t_ant - S_t_ant_new2; 
 % plot(abs(S_t_ant(:)));
 % set(gca, 'YScale', "Linear");
 % yline(TH_abs(1));
 % hold on;
 % plot(abs(test(:)));


% freq_response = abs(fft(fir_coeff_sinc, Nfft));  % Модуль спектра фильтра
% freq_axis = linspace(0, 1, Nfft);  % Частотная ось (нормализованная от 0 до 1)


% figure;
% stem(fir_coeff_sinc, 'filled');
% title('во временной области');
% xlabel('Время');
% ylabel('Коэффициенты фильтра');
% grid on;
% 
% % Построение амплитудного спектра
% figure;
% plot(freq_axis, freq_response);
% title('');
% xlabel('Нормализованная частота');
% ylabel('Амплитуда');
% grid on;
% end


