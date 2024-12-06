function [Xtransformed] = ft_transform(X,mode,scenario)
% X [Nsym, Nsc, Ntx]

Ntx = scenario.Ntx;
Nfft = scenario.Nfft;
Nsym = scenario.Nsym;
Nsc = scenario.Nsc;

shift = Nsc/2;

if mode == "f2t"
    X_f = zeros(Nsym,Nfft,Ntx);
    X_f(:,1:Nsc,:) = X;
    X_f = circshift(X_f, -shift, 2);
    
    Xtransformed = ifft(X_f,Nfft,2) * sqrt(Nfft);
elseif mode == "t2f"
    X_f = fft(X,Nfft,2) / sqrt(Nfft); % TODO rename this
    X_f = circshift(X_f, shift,2);
    Xtransformed = X_f(:,1:Nsc,:);
end