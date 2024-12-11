% function [dX] = SINC_paper(X_t, A, N_used, TH, max_EVM, Nfft, scen)
% Nsym=size(X_t,1); % Number of OFDM symbols in the initial time-domain signal
% 
% Ndac=scen.Ndac; % Number of digital-to-analog converters
% Nant=scen.Ntx; % Number of antennas (> Ndac, HBF architecture)
% 
% N_samples=16; % Number of points at the interval 
%               % (1 interval - 1 peak per iteration)
% 
% coef1=0; % Use only peaks in inversion
% coef2=Nant/Ndac; % Use peaks and zeros in inversion
%                  % Default ratio for no power loss for coef2: Nant/Ndac 
%                  % (Since a peak is a delta function we save only Ndac/Nant
%                  % power when cut DAC domain Nant/Ndac times)
% 
% Nzero=Nfft-N_used; % padding for IFFT calculation
% 
% S_t_ant = permute(X_t, [1 3 2]); % Initial signal to reduce PAPR in
%                                  % Size: [Nsym, Nant, Nfft]
% 
% 
% % signal power and root mean square
% mean_power=mean( abs(S_t_ant(:,:,:)).^2,'all');
% 
% % Tunable threasholds (step-by-step reduction can be implemented)
% % with increasing the TH, the number of selected peaks is reduced,
% % resulting in less PAPR reduction, but less Error vector magnitude (EVM)
% %TH = [2.2 1.3];  % SET the THREASHOLD VALUE !!!
% TH_abs=TH*sqrt( mean_power );
% N_iter = length(TH);
% 
% 
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
% 
% S_t_canc=zeros(Nsym,Nant,Nfft);
% S_t_ant_new=zeros(Nsym,Nant,Nfft);
% S_t_ant_canc_peak=zeros(Nsym,Nant,Nfft);
% 
% % Number of iterations for peaks positions search per interval
% M_iter = 2;
% 
% %%%%%%%%%%%%%%New way to find intervals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %find peaks to reduce
% for i1=1:Nsym
%     % fprintf("i1 =%d\n", i1);
%     for i2=1:Nant
%         % fprintf("i2 =%d\n", i2);
%         % Extract signal at particular symbol at particular antenna
%         S_t = squeeze(S_t_ant(i1,i2,:)).';
%         % Извлекаем интервалы для текущего символа и антенны
%         min_inds = find_intervals(S_t);
% 
%         % Пошаговое подавление пиков
%         for j = 1:N_iter
%             S_t_canc_tmp = zeros(1, Nfft);
% 
%             % Проходим по всем интервалам
%             for k = 1:length(min_inds) - 1
%                 interval_start = min_inds(k);
%                 interval_end = min_inds(k + 1) - 1; % Интервалы закрыты
% 
%                 % Извлекаем сигнал на интервале
%                 signal = S_t(interval_start:interval_end);
% 
%                 % Находим максимум и его индекс
%                 [Max_value, Indx] = max(abs(signal));
% 
%                 % Проверяем порог
%                 if Max_value > TH_abs(j)
%                     % Амплитуда подавления
%                     sinc_Ampl = signal(Indx) * (1 - TH_abs(j) / Max_value);
% 
%                     % Глобальная позиция пика
%                     sinc_shift = interval_start + Indx - 1;
% 
%                     % Суммируем синусоиды для подавления
%                     S_t_canc_tmp = S_t_canc_tmp + sinc_Ampl * circshift(SINC_t, [0, sinc_shift - 1]);
% 
%                     % Запоминаем подавляемый пик
%                     S_t_ant_canc_peak(i1, i2, sinc_shift) = ...
%                         S_t_ant_canc_peak(i1, i2, sinc_shift) + sinc_Ampl;
%                 end
%             end
% 
%             % Обновляем сигнал после подавления
%             S_t = S_t - S_t_canc_tmp;
%             S_t_canc(i1, i2, :) = squeeze(S_t_canc(i1, i2, :)) + S_t_canc_tmp.';
%         end
% 
%         % Сохраняем модифицированный сигнал
%         S_t_ant_new(i1, i2, :) = S_t;
% 
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Old variant 
% 
% % for i1=1:Nsym
% % 
% %     for i2=1:Nant
% % 
% %         % Extract signal at particular symbol at particular antenna
% %         S_t = squeeze(S_t_ant(i1,i2,:)).';
% % 
% %         % PAPR reduction algorithm
% %         for j=1:N_iter
% %             S_t_canc_tmp=zeros(1,Nfft);
% %             for p = 1:M_iter
% %                 Nintervals = Nfft/N_samples; 
% %                  % Создаём новый рисунок
% %                 % figure;
% %                 % plot(abs(S_t), 'b-', 'LineWidth', 1.5); % Амплитуда сигнала
% %                 % hold on;
% %                 % 
% %                 % % Отображение интервалов
% %                 % for k = 1:Nintervals
% %                 %     % Границы текущего интервала
% %                 %     interval_start = (k - 1) * N_samples + 1;
% %                 %     % interval_end = min(k * N_samples, Nfft); % Учитываем возможный выход за пределы Nfft
% %                 % 
% %                 %     % Рисуем границы интервалов
% %                 %     xline(interval_start, 'r--', 'LineWidth', 1.2); % Красная линия — начало интервала
% %                 %     % xline(interval_end, 'g--', 'LineWidth', 1.2);   % Зелёная линия — конец интервала
% %                 % end
% %                 % 
% %                 % % Оформление графика
% %                 % title(sprintf('Original intervals', i1, i2));
% %                 % xlabel('Index');
% %                 % ylabel('Amplitude');
% %                 % legend('Signal', 'Bound');
% %                 % hold off;
% %                 % Peaks positions and amplidute per interval
% %                 sinc_Ampl = zeros(Nintervals,1);
% %                 sinc_shift = zeros(Nintervals,1);
% %                 for k=1:Nintervals
% % 
% %                     % Extract signal at particular interval
% %                     signal = S_t((k-1)*N_samples+1:k*N_samples);
% % 
% %                     [Max_value,Indx] = max(abs(signal));
% % 
% %                     % Check for threshold exceeding and peak parameters set
% %                     if Max_value>TH_abs(j)
% % 
% %                         % We are interested in peak value above the threshold
% %                         sinc_Ampl(k) = signal(Indx)*(1 - TH_abs(j)/Max_value);
% % 
% %                         % Global (among Nfft points) peak position
% %                         sinc_shift(k) = Indx+(k-1)*N_samples;
% %                     else
% %                         sinc_Ampl(k)=0;
% %                         sinc_shift(k)=1;
% %                     end
% % 
% %                     % Sum of sincs canceling found peaks
% %                     S_t_canc_tmp = S_t_canc_tmp + sinc_Ampl(k)*circshift(SINC_t,[0 sinc_shift(k)-1]);
% % 
% %                     % Sum of delta functions canceling found peaks 
% %                     % (they are transformed to DAC domain)
% %                     S_t_ant_canc_peak(i1,i2,Indx+(k-1)*N_samples) = S_t_ant_canc_peak(i1,i2,Indx+(k-1)*N_samples) + sinc_Ampl(k); 
% %                 end
% % 
% %                 S_t = S_t - S_t_canc_tmp;
% %                 S_t_canc(i1,i2,:) = squeeze(S_t_canc(i1,i2,:)) + S_t_canc_tmp.';
% % 
% %             end
% %         end
% % 
% %         % Modified signal after per-antenna PAPR reduction
% %         S_t_ant_new(i1,i2,:) = S_t;
% % 
% %     end
% % end
% 
% % % Now we need to check what happens in DAC domain
% S_t_dac_canc_sig=zeros(Nsym,Ndac,Nfft);
% S_t_dac_canc_peak = zeros(Nsym,Ndac,Nfft);
% S_t_ant_canc_peak_new=zeros(Nsym,Nant,Nfft);
% 
% for i1=1:Nsym
%     for i3=1:Nfft
%         % Extract signal at all antennas
%         sig = squeeze(S_t_ant_canc_peak(i1,:,i3));
% 
%         % Check for found peaks
%         indx = find(abs(sig));
% 
%         if sum(indx)>0
%             A_new = A(indx,:);
%             sig_new = sig(indx);
%             s1 = pinv(A_new) * sig_new.'; % peaks only transform
%             s2 = pinv(A) * sig.'; % peaks + zeros transform
%             S_t_dac_canc_peak(i1,:,i3)=coef1*s1+coef2*s2;
%         else
%             S_t_dac_canc_peak(i1,:,i3) = zeros(1,Ndac);
%         end
% 
%     end
%     % Go back to antenna domain 
%     S_t_ant_canc_peak_new(i1,:,:) = A * squeeze(S_t_dac_canc_peak(i1,:,:));
% 
%     % Convolve peaks with sinc functions in DAC domain to satisfy spectrum
%     % mask
%     for i2=1:Ndac
%         s_t = squeeze(S_t_dac_canc_peak(i1,i2,:));
% 
%         S_t_dac_canc_sig(i1,i2,:) = s_t.' * SINC_mtx;
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check for EVM and reduce CFR noise proportionally
% EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2,'all') / sum(abs(S_t_ant).^2,'all'));
% S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx/max_EVM,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Generate the final signal at antennas to analyze PAPR
% S_t_ant_new2 = zeros(Nsym,Nant,Nfft);
% for i1=1:Nsym
%     S_t_ant_new2(i1,:,:) = A * squeeze(S_t_dac_canc_sig(i1,:,:));
% end
% 
% dS = S_t_ant_new2;
% dX = permute(dS, [1,3,2]);

% 
function [dX] = SINC_paper(X_t, A, N_used, TH, max_EVM, Nfft, scen)
    Nsym = size(X_t, 1); % Количество символов OFDM в исходном временном сигнале
    Ndac = scen.Ndac; % Количество цифро-аналоговых преобразователей
    Nant = scen.Ntx; % Количество антенн (> Ndac, HBF архитектура)

    N_samples = 16; % Количество точек на интервал (1 интервал - 1 пик на итерацию)
    coef1 = 0; % Использовать только пики в инверсии
    coef2 = Nant / Ndac; % Использовать пики и нули в инверсии
    Nzero = Nfft - N_used; % Паддинг для расчета IFFT

    S_t_ant = permute(X_t, [1 3 2]); % Исходный сигнал для уменьшения PAPR (Размер: [Nsym, Nant, Nfft])

    % Мощность сигнала и среднеквадратичное значение
    mean_power = mean(abs(S_t_ant(:,:,:)).^2, 'all');

    % Пороговые значения
    TH_abs = TH * sqrt(mean_power);
    N_iter = length(TH);

    % Параметры фильтра Баттерворта
    n = 4; % Порядок фильтра
    Fc = 0.2; % Частота среза (нормализованная частота от 0 до 1)
    [b, a] = butter(n, Fc, 'low'); % Фильтр Баттерворта низких частот

    % Генерация SINC сигнала
    SINC_f = circshift([ones(1, N_used), zeros(1, Nzero)], [0, -(N_used) / 2]);
    SINC_t = ifft(SINC_f) * sqrt(Nfft);
    SINC_t = SINC_t / SINC_t(1); % Нормализация амплитуды для совпадения с пиковой амплитудой

    % Генерация матрицы SINC для различных временных позиций
    SINC_mtx = zeros(Nfft, Nfft);
    for j = 1:Nfft
        SINC_mtx(j, :) = circshift(SINC_t, [0, j - 1]);
    end

    % Инициализация переменных
    S_t_canc = zeros(Nsym, Nant, Nfft);
    S_t_ant_new = zeros(Nsym, Nant, Nfft);
    S_t_ant_canc_peak = zeros(Nsym, Nant, Nfft);

    % Количество итераций для поиска пиков в интервале
    M_iter = 2;

    % Найдем пики для подавления
    for i1 = 1:Nsym
        for i2 = 1:Nant
            % Извлекаем сигнал для текущего символа и антенны
            S_t = squeeze(S_t_ant(i1, i2, :)).';

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

                        % Применяем фильтр Баттерворта (фильтруем сигнал для подавления пиков)
                        filtered_signal = filter(b, a, circshift(SINC_t, [0, sinc_shift - 1]));

                        % Добавляем откорректированный сигнал
                        S_t_canc_tmp = S_t_canc_tmp + sinc_Ampl * filtered_signal;

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

    % Теперь проверим, что происходит в DAC домене
    S_t_dac_canc_sig = zeros(Nsym, Ndac, Nfft);
    S_t_dac_canc_peak = zeros(Nsym, Ndac, Nfft);
    S_t_ant_canc_peak_new = zeros(Nsym, Nant, Nfft);

    for i1 = 1:Nsym
        for i3 = 1:Nfft
            % Извлекаем сигнал для всех антенн
            sig = squeeze(S_t_ant_canc_peak(i1, :, i3));

            % Проверяем найденные пики
            indx = find(abs(sig));

            if sum(indx) > 0
                A_new = A(indx, :);
                sig_new = sig(indx);
                s1 = pinv(A_new) * sig_new.'; % Только пики в преобразовании
                s2 = pinv(A) * sig.'; % Пики + нули в преобразовании
                S_t_dac_canc_peak(i1, :, i3) = coef1 * s1 + coef2 * s2;
            else
                S_t_dac_canc_peak(i1, :, i3) = zeros(1, Ndac);
            end
        end
        % Возвращаемся в антенный домен
        S_t_ant_canc_peak_new(i1, :, :) = A * squeeze(S_t_dac_canc_peak(i1, :, :));

        % Применяем пики с фильтром Баттерворта в DAC домене для удовлетворения спектральной маски
        for i2 = 1:Ndac
            s_t = squeeze(S_t_dac_canc_peak(i1, i2, :));
            S_t_dac_canc_sig(i1, i2, :) = s_t.' * SINC_mtx; % Применяем SINC для моделирования спектра
        end
    end

    % Проверка на EVM и пропорциональное уменьшение шума CFR
    EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2, 'all') / sum(abs(S_t_ant).^2, 'all'));
    S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx / max_EVM, 1);

    % Генерируем финальный сигнал для антенн для анализа PAPR
    S_t_ant_new2 = zeros(Nsym, Nant, Nfft);
    for i1 = 1:Nsym
        S_t_ant_new2(i1, :, :) = A * squeeze(S_t_dac_canc_sig(i1, :, :));
    end

    % Результирующий сигнал
    dS = S_t_ant_new2;
    dX = permute(dS, [1, 3, 2]);
end
