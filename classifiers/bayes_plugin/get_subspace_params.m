function P = get_subspace_params(adjacency_matrices,constants,subspace)

% compute parameters for choosing subspaces
if isfield(subspace,'naive_bayes') || isfield(subspace,'incoherent') || isfield(subspace,'incoherent')
    P.d_inc=abs(P.E0-P.E1);
end

if isfield(subspace,'cor_incoherent')
    P.d_zinc = abs(P.E0./sqrt(P.E0.*(1-P.E0)) - P.E1./sqrt(P.E1.*(1-P.E1)));
end

if isfield(subspace,'fisher_incoherent')
    P.d_pval = run_get_fisher_pvals(adjacency_matrices,constants);
end
