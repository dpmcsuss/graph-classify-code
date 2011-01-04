function [Lhat params subspace Lsem] = run_naive_bayes(Atrn,Gtrn,Atst,Gtst,subspace)
% get mle and do naive bayes classifier


params = get_params_mle(Atrn,Gtrn);

if isfield(subspace,'fisher_incoherent')
    if length(subspace.fisher_incoherent)<2
        params.d_pval = run_get_fisher_pvals(As,constants);
    end
end

subspace = get_subspace(params,subspace);

[Lhat Lsem] = naive_bayes_classify(Atst,Gtst.ys,params,subspace);

end