function cuts = fcp_potential_breaks(cfg)
cuts=[];
if isfield(cfg,'q_breaks'), cuts=cfg.q_breaks(:); end
end
