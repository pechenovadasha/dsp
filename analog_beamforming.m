function Fa = analog_beamforming(H, scenario)
Ntx = scenario.Ntx;
Ndac = scenario.Ndac;

% Same as beamspace selection (choose powerful directions from FFT mtx)
F = dftmtx(Ntx);

H = reshape(H, [], Ntx);

P = H * F;

P = sum(abs(P).^2 , 1);

[~, idx] = sort(P, 'descend');

Fa = F(:,idx(1:Ndac)) / sqrt(Ntx);