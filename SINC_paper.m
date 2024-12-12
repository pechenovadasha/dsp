function [dX] = SINC_paper(X_t, A, N_used, TH, max_EVM, Nfft, scen)
        
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
    mean_power = mean(abs(S_t_ant(:,:,:)).^2, 'all');


    TH_abs = TH * sqrt(mean_power);
    N_iter = length(TH);

    % Параметры фильтра Баттерворта
    n = 4; % Порядок фильтра
    Fc = 0.2; % Частота среза (нормализованная частота от 0 до 1)
    [b, a] = butter(n, Fc, 'low'); % Фильтр Баттерворта низких частот


    SINC_f = circshift([ones(1, N_used), zeros(1, Nzero)], [0, -(N_used) / 2]);
    SINC_t = ifft(SINC_f) * sqrt(Nfft);
    SINC_t = SINC_t / SINC_t(1); 

    
    SINC_mtx = zeros(Nfft, Nfft);
    for j = 1:Nfft
        SINC_mtx(j, :) = circshift(SINC_t, [0, j - 1]);
    end


    S_t_canc = zeros(Nsym, Nant, Nfft);
    S_t_ant_new = zeros(Nsym, Nant, Nfft);
    S_t_ant_canc_peak = zeros(Nsym, Nant, Nfft);

    for i1 = 1:Nsym
        for i2 = 1:Nant
  
            S_t = squeeze(S_t_ant(i1, i2, :)).';

            min_inds = find_intervals(S_t);

        
            for j = 1:N_iter
                S_t_canc_tmp = zeros(1, Nfft);

       
                for k = 1:length(min_inds) - 1
                    interval_start = min_inds(k);
                    interval_end = min_inds(k + 1) - 1; 

                    signal = S_t(interval_start:interval_end);

             
                    [Max_value, Indx] = max(abs(signal));

                   
                    if Max_value > TH_abs(j)
                     
                        sinc_Ampl = signal(Indx) * (1 - TH_abs(j) / Max_value);

                        
                        sinc_shift = interval_start + Indx - 1;

                      
                        filtered_signal = filter(b, a, circshift(SINC_t, [0, sinc_shift - 1]));

                       
                        S_t_canc_tmp = S_t_canc_tmp + sinc_Ampl * filtered_signal;

                        
                        S_t_ant_canc_peak(i1, i2, sinc_shift) = ...
                            S_t_ant_canc_peak(i1, i2, sinc_shift) + sinc_Ampl;
                    end
                end

               
                S_t = S_t - S_t_canc_tmp;
                S_t_canc(i1, i2, :) = squeeze(S_t_canc(i1, i2, :)) + S_t_canc_tmp.';
            end

            
            S_t_ant_new(i1, i2, :) = S_t;

        end
    end


    S_t_dac_canc_sig = zeros(Nsym, Ndac, Nfft);
    S_t_dac_canc_peak = zeros(Nsym, Ndac, Nfft);
    S_t_ant_canc_peak_new = zeros(Nsym, Nant, Nfft);

    for i1 = 1:Nsym
        for i3 = 1:Nfft
            
            sig = squeeze(S_t_ant_canc_peak(i1, :, i3));

            
            indx = find(abs(sig));

            if sum(indx) > 0
                A_new = A(indx, :);
                sig_new = sig(indx);
                s1 = pinv(A_new) * sig_new.';
                s2 = pinv(A) * sig.';
                S_t_dac_canc_peak(i1, :, i3) = coef1 * s1 + coef2 * s2;
            else
                S_t_dac_canc_peak(i1, :, i3) = zeros(1, Ndac);
            end
        end
       
        S_t_ant_canc_peak_new(i1, :, :) = A * squeeze(S_t_dac_canc_peak(i1, :, :));

  
        for i2 = 1:Ndac
            s_t = squeeze(S_t_dac_canc_peak(i1, i2, :));
            S_t_dac_canc_sig(i1, i2, :) = s_t.' * SINC_mtx; 
        end
    end

    EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2, 'all') / sum(abs(S_t_ant).^2, 'all'));
    S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx / max_EVM, 1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Check for EVM and reduce CFR noise proportionally
% EVM_approx = sqrt(sum(abs(S_t_dac_canc_sig).^2,'all') / sum(abs(S_t_ant).^2,'all'));
% S_t_dac_canc_sig = S_t_dac_canc_sig / max(EVM_approx/max_EVM,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Generate the final signal at antennas to analyze PAPR

    S_t_ant_new2 = zeros(Nsym, Nant, Nfft);
    for i1 = 1:Nsym
        S_t_ant_new2(i1, :, :) = A * squeeze(S_t_dac_canc_sig(i1, :, :));
    end

    dS = S_t_ant_new2;
    dX = permute(dS, [1, 3, 2]);
end
