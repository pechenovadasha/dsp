function X = digital_twin(S,Fa,Fd,scenario)

Ntx = scenario.Ntx;
Nsc = scenario.Nsc;
Nsym = scenario.Nsym;

X = zeros(Nsym,Nsc,Ntx);

for sym = 1:Nsym
    for sc = 1:Nsc
        X(sym,sc,:) = Fa*squeeze(Fd(sc,:,:))*squeeze(S(sym,sc,:));
    end
end
