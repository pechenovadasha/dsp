function [Xm] = PAPR_reduction(X,Fa,scenario,ML_coef)

max_evm = scenario.max_evm;

Ntx = scenario.Ntx;
Nsc = scenario.Nsc;
Nfft = scenario.Nfft;
Nsym = scenario.Nsym;

X_t = ft_transform(X,"f2t",scenario);

threshold = db2mag(ML_coef); % dB
algo = scenario.PAPR_algo;
switch algo
    case 1 % For original signal calculation
        dX_t = zeros(Nsym, Nfft, Ntx); 
    case 2
        [dX_t] = SINC_2D(X_t, Fa, Nsc, threshold, max_evm, Nfft,scenario);
    case 3

        [dX_t] = SINC_paper(X_t, Fa, Nsc, threshold, max_evm, Nfft, scenario);
end
dX = ft_transform(dX_t,"t2f",scenario);

Xm = zeros(Nsym,Nsc,Ntx);
for sym = 1:Nsym
    dZ = squeeze(dX(sym,:,:))*pinv(Fa).'; % Put it into STR to avoid this code
    dX_feasible = dZ*Fa.';

    Xm(sym,:,:) = squeeze(X(sym,:,:)) - dX_feasible;
end

