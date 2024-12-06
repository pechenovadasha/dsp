function metrics = tester(scenario, ML_coef)
Nsc = scenario.Nsc;
Nue = scenario.Nue;
Ntx = scenario.Ntx;

% Generate data for transmission
S = generate_data(scenario);

% Load channel (ideal CE assumed)
H = (randn(Nsc,Nue,Ntx)+1i*randn(Nsc,Nue,Ntx))/sqrt(2);

% Compute analog beamforming (beamspace selection)
Fa = analog_beamforming(H,scenario);

% Compute digital beamforming (no precoding - eye)
Fd = repmat(reshape(eye(Nue),[1,Nue,Nue]),[Nsc,1,1]);

% Compute TX signals before PAPR reduction
[X] = digital_twin(S,Fa,Fd,scenario);

% PAPR compensated signal
[Xm] = PAPR_reduction(X,Fa,scenario, ML_coef);

% Compute PAPR, CCDF
metrics = compute_metrics(Xm, scenario);

