function [alpha,alpha_w,alpha_q,info] = fcp_parameter_rule(cfg,model,relative)
scale=max(model.singular_values(1)^2,realmin);
candidates=[cfg.alpha_factor*relative^cfg.alpha_power,cfg.alpha_min,100*eps*scale];
[alpha,active]=max(candidates);
alpha_w_rule=cfg.alpha_w_factor*relative^cfg.alpha_w_power;
alpha_w=max(alpha_w_rule,cfg.alpha_w_min);
alpha_q_rule=cfg.alpha_q_factor*relative^cfg.alpha_q_power;
alpha_q=max([alpha_q_rule,cfg.alpha_q_min,100*eps]);
info.alpha_lower_bound_active=active>1;
info.alpha_w_lower_bound_active=alpha_w>alpha_w_rule;
info.alpha_q_lower_bound_active=alpha_q>alpha_q_rule;
end
