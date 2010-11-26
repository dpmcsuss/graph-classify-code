function [Lhat subspace] = run_ind_edge(Atrn,Gtrn,Atst,Gtst,subspace)
% get mle and do naive bayes classifier


params = get_params_mle(Atrn,Gtrn);

if isfield(subspace,'incoherent')
    if isempty(subspace.incoherent)
        subspace.incoherent=get_inc_edges(params.d_pos,Gtrn.n);
    end
end
        
if isfield(subspace,'cor_incoherent')
    if isempty(subspace.cor_incoherent)
        subspace.cor_incoherent=get_inc_edges(params.d_opt,Gtrn.n);
    end
end

if isfield(subspace,'coherent')
    if isempty(subspace.coherent)
        subspace.coherent=get_max_edges(params.d_pos,7);
    end
end


Lhat = ind_edge_classify(Atst,Gtst.ys,params,subspace);

end