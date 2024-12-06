function X = generate_data(scenario)

Nsym = scenario.Nsym;
Nue = scenario.Nue;
Nsc = scenario.Nsc;
QAM_order = scenario.QAM_order;
QAM_points = 2^QAM_order;

data_bits = randi([0 1], [QAM_order, Nsym*Nsc*Nue]);

data_syms = bit2int_(data_bits,QAM_order);
data_syms = qammod(data_syms,QAM_points,'UnitAveragePower', true); % Gray coding, phase offset = 0;

X = reshape(data_syms, [Nsym, Nsc, Nue]);