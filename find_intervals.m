function intervals = find_intervals(x_t)
    % Получаем размеры входного сигнала
   [Nsym, Nfft, Nant] = size(x_t);

    % Массив для хранения интервалов

    intervals = [];
    % Проходим по всем символам и антеннам
    for i1 = 1:Nsym
        for i2 = 1:Nant
            % Извлекаем сигнал для текущего символа и антенны
            signal = squeeze(x_t(i1, :, i2)); 

            % Работаем с амплитудой сигнала
            abs_signal = abs(signal);

            % Находим локальные минимумы вручную
            min_inds = [];
            for k = 2:Nfft-1
                if abs_signal(k) < abs_signal(k-1) && abs_signal(k) < abs_signal(k+1)
                    min_inds = [min_inds, k]; % Добавляем индекс минимума
                end
            end

            % Убедимся, что у нас есть хотя бы один минимум
            if isempty(min_inds)
                min_inds = [1, Nfft]; % Если минимумов нет, берем весь сигнал
            else
                % Добавляем границы сигнала, если их нет в списке минимумов
                if min_inds(1) ~= 1
                    min_inds = [1, min_inds];
                end
                if min_inds(end) ~= Nfft
                    min_inds = [min_inds, Nfft];
                end
            end

          
            % % Разбиение на интервалы, где границы - точки минимума
            num_intervals = length(min_inds) - 1;
            intervals_for_current = cell(num_intervals, 1);

            % % Создаем новый рисунок для отображения
            % figure;
            % plot(abs_signal, 'b-', 'LineWidth', 1.5); % Сигнал в виде линии
            % hold on;
            % 
            % for j = 1:num_intervals
            %     interval_start = min_inds(j);
            %     interval_end = min_inds(j + 1);
            % 
            %     % Сохраняем текущий интервал
            %     intervals_for_current{j} = signal(interval_start:interval_end);
            % 
            %     % Отображаем границы интервалов на графике
            %     xline(interval_start, 'r--', 'LineWidth', 1.2); % Красная линия - начало интервала
            %     xline(interval_end, 'g--', 'LineWidth', 1.2);   % Зеленая линия - конец интервала
            % end
            % 
            % % Добавляем обозначения
            % title('Division by minimums');
            % xlabel('Index');
            % ylabel('Amplitude');
            % legend('Signal', 'Start of interval', 'End of interval');
            % hold off;
        end
    end
    intervals = min_inds;
end
