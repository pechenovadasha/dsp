function [CCDF] = PAPR_analyzer(S_t,scenario)
    
    Nsym = scenario.Nsym;
    Nfft = scenario.Nfft;
    
    % If input is in frequency domain
    if size(S_t,2) < Nfft
        S_t = ft_transform(S_t,"f2t",scenario);
    end

    
    PAPR_dB = zeros(Nsym,1);
    PAPR_full = [];
    % Compute PAPR per symbol + aggregate all samples to array for CCDF
    
    P = abs(S_t).^2;
    mean_P = mean(P,'all');
    PAPR_full = P./mean_P;
    PAPR_full = pow2db(PAPR_full(:));
    
    papr = scenario.papr_ccdf;
    for i = 1:length(papr)
        y(i) = sum(PAPR_full >= papr(i))/length(PAPR_full); % # of samples exceeding papr(i)
    end
 	CCDF = y;
end