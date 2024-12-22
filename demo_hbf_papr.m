rng(42)

scenario.Ntx = 64; % transmit antennas
scenario.Ndac = 16; % digital-to-analog converters
scenario.Nue = 16; % single-antenna users
scenario.Nsc = 192;

scenario.Nsym = 2; % symbols in time
scenario.Nfft = 2048;

scenario.QAM_order = 8;
scenario.max_evm = 0.125;


% PAPR_thresholds = [8 8 8]; % peak threshold in dB per iteration
PAPR_thresholds = [8 8 8];
scenario.papr_ccdf = 5:0.1:15;

scenario.N_scen = 12; % Number of scenarios to simulate 
                     % (ccdf 10^-6 should be achieved)

scenario.PAPR_algo = 1;
metrics = scenarios_validation(scenario,PAPR_thresholds);
semilogy(scenario.papr_ccdf, mean(reshape([metrics.CCDF],numel(scenario.papr_ccdf),scenario.N_scen),2));
hold on

scenario.PAPR_algo = 3;
metrics = scenarios_validation(scenario,PAPR_thresholds);
semilogy(scenario.papr_ccdf, mean(reshape([metrics.CCDF],numel(scenario.papr_ccdf),scenario.N_scen),2));

legend("Original signal", "PAPR reduction")


function metrics = scenarios_validation(scenario, ML_coef)
    
N_scen = scenario.N_scen;
metrics = struct('CCDF', cell(N_scen, 1));

for index_i = 1:N_scen 
    rng(index_i);

    scen = scenario;

    scen.index = index_i + double(index_i>23); % 24 is very bad scenario, smth wrong with precoding?
    
    metrics(index_i) = tester(scen, ML_coef);
end

end
