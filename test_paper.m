function [dX] = test_paper(X_t, A, N_used, TH, max_EVM, Nfft, scen)
    Nsym = size(X_t, 1);  % Количество символов OFDM в исходном временном сигнале
    Ndac = scen.Ndac;     % Количество цифро-аналоговых преобразователей
    Nant = scen.Ntx;      % Количество антенн (> Ndac, HBF архитектура)

    Nzero = Nfft - N_used;  % Дополнение нулями для расчёта IFFT

    S_t_ant = permute(X_t, [1 3 2]);  % Исходный сигнал для снижения PAPR
                                      % Размер: [Nsym, Nant, Nfft]

    % Средняя мощность сигнала и среднеквадратическое значение
    mean_power = mean(abs(S_t_ant(:,:,:)).^2, 'all');

    % Устанавливаем пороги на основе среднего уровня мощности
    TH_abs = TH * sqrt(mean_power);
    N_iter = length(TH);

    % Генерация SINC-функции
    SINC_f = circshift([ones(1, N_used) zeros(1, Nzero)], [0 -(N_used)/2]);
    SINC_t = ifft(SINC_f) * sqrt(Nfft);
    SINC_t = SINC_t / SINC_t(1);  % Нормализация амплитуды

    % Создание матрицы сдвигов SINC-функции
    SINC_mtx = zeros(Nfft, Nfft);
    for j = 1:Nfft
        SINC_mtx(j, :) = circshift(SINC_t, [0 j-1]);
    end

    S_t_canc = zeros(Nsym, Nant, Nfft);

    % Основной процесс подавления пиков
    for i1 = 1:Nsym
        for i2 = 1:Nant
            % Извлекаем сигнал для текущего символа и антенны
            S_t = squeeze(S_t_ant(i1, i2, :)).';

            % Находим интервалы для подавления пиков
            min_inds = find_intervals(S_t);

            % Оптимизация для каждого порога
            for j = 1:N_iter
                % Формулировка задачи оптимизации в CVX
                cvx_begin quiet
                    variable S_m(Ndac, Nfft) complex;  % Искомый сигнал в сжатом пространстве
                    minimize( norm(S_t - A * S_m * SINC_mtx, 'fro') )
                    subject to
                        abs(S_t - A * S_m * SINC_mtx) < TH_abs(j);  % Пороговое ограничение
                cvx_end

                % Обновляем сигнал на текущей итерации
                S_t = S_t - A * S_m * SINC_mtx;
                S_t_canc(i1, i2, :) = squeeze(S_t_canc(i1, i2, :)) + (A * S_m * SINC_mtx).';
            end
        end
    end

    % Проверка EVM и нормализация сигнала
    S_t_dac_canc_sig = zeros(Nsym, Ndac, Nfft);
    for i1 = 1:Nsym
        for i2 = 1:Ndac
            s_t = squeeze(S_t_canc(i1, i2, :));
            S_t_dac_canc_sig(i1, i2, :) = s_t.' * SINC_mtx;
        end
    end

    % Оценка EVM
    EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2, 'all') / sum(abs(S_t_ant).^2, 'all'));
    S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx / max_EVM, 1);

    % Генерация итогового сигнала на антеннах
    S_t_ant_new2 = zeros(Nsym, Nant, Nfft);
    for i1 = 1:Nsym
        S_t_ant_new2(i1, :, :) = A * squeeze(S_t_dac_canc_sig(i1, :, :));
    end

    % Вывод результата
    dS = S_t_ant_new2;
    dX = permute(dS, [1, 3, 2]);
end
