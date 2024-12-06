function metrics = compute_metrics(Xm, scenario)

[metrics.CCDF] = PAPR_analyzer(Xm,scenario);

