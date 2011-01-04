function [subspace params] = get_subspace(params,subspace)
% this function gets subspaces for whatever subspace fields are specified,
% assuming that they are one of the subspaces for which we have built in
% algorithms


n=length(params.E0);

% get subspaces
if isfield(subspace,'naive_bayes') % use all edges for which there is a class-conditional difference
    subspace.naive_bayes = find(params.d_pos>0);
end

if isfield(subspace,'incoherent')
    if isempty(subspace.incoherent)
        subspace.incoherent=get_inc_edges(params.d_pos,n);
    elseif numel(subspace.incoherent)==1
        subspace.incoherent=get_inc_edges(params.d_pos,subspace.incoherent);
    end
end

if isfield(subspace,'cor_incoherent')   % incoherent edges ranked by "z-score"
    if isempty(subspace.cor_incoherent)
        subspace.cor_incoherent=get_inc_edges(params.d_opt,n);
    elseif numel(subspace.cor_incoherent)==1
        subspace.cor_incoherent=get_inc_edges(params.d_opt,subspace.cor_incoherent);
    end
end

if isfield(subspace,'fisher_incoherent')
    if isempty(subspace.fisher_incoherent)
        subspace.fisher_incoherent=get_inc_edges(1-params.d_pval,n);
    elseif numel(subspace.fisher_incoherent)==1
        subspace.fisher_incoherent=get_inc_edges(1-params.d_pval,subspace.fisher_incoherent);
    end
end

if isfield(subspace,'coherent')
    if isempty(subspace.coherent)
        subspace.coherent=get_max_edges(params.d_pos,round(sqrt(n)));
    elseif numel(subspace.coherent)==1
        subspace.coherent=get_max_edges(params.d_pos,subspace.coherent);
    end
end
