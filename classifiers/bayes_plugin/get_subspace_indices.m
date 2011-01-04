function [subspace params] = get_subspace_indices(params,subspace,As,constants)
% this function gets subspaces for whatever subspace fields are specified,
% assuming that they are one of the subspaces for which we have built in
% algorithms

for i=1:length(subspace)                % for each subspace specified

    % ensure 'indices' field is present
    if ~isfield(subspace(i),'indices'), subspace(i).indices=[]; end
    
    % get indices if they are not specified
    if isempty(subspace(i).indices)     % if indices are not pre-specified
        switch subspace(i).name         % get indices using supplied method
            
            case 'naive_bayes'          % use all edges that have class conditional different distributions
                if ~isfield(params,'d_inc')
                    params.d_inc=abs(params.E0-params.E1);
                end
                subspace(i).indices = find(params.d_inc>1e-4);
                subspace(i).size    = length(subspace(i).indices);

            case 'incoherent'           % use absolute difference to rank edges
                if ~isfield(params,'d_inc')
                    params.d_inc=abs(params.E0-params.E1);
                end
                subspace(i).size    = min(length(find(params.d_inc>0)),subspace(i).size);
                subspace(i).indices = get_inc_edges(params.d_inc,subspace(i).size);
                
            case 'incoherent_z'         % use z-score to rank edges
                if ~isfield(params,'d_zinc')
                    params.d_zinc = abs(params.E0./sqrt(params.E0.*(1-params.E0)) - params.E1./sqrt(params.E1.*(1-params.E1)));
                end
                subspace(i).size    = min(length(find(params.d_zinc>0)),subspace(i).size);
                subspace(i).indices = get_inc_edges(params.d_zinc,subspace(i).size);
                
            case 'incoherent_f'         % use exact fisher's text to rank edges
                if ~isfield(params,'d_pval')
                    params.d_pval = run_get_fisher_pvals(As,constants);
                end
                subspace(i).size    = min(length(find(params.d_pval<1)),subspace(i).size);
                subspace(i).indices = get_inc_edges(1-params.d_pval,subspace(i).size);
                
            case 'coherent'             % use degree to rank vertices
                if ~isfield(params,'d_inc')
                    params.inc=abs(params.E0-params.E1);
                end
                subspace(i).indices = get_max_edges(params.d_inc,subspace(i).size);
                
            case 'coherent_z'           % use summed z-scores to rank vertices
                if ~isfield(params,'d_zinc')
                    params.d_zinc = abs(params.E0./sqrt(params.E0.*(1-params.E0)) - params.E1./sqrt(params.E1.*(1-params.E1)));
                end
                subspace(i).indices = get_max_edges(params.z_inc,subspace(i).size);
                
            case 'coherent_f'           % use 1-summed pvals to rank vertices
                if ~isfield(params,'d_pval')
                    params.d_pval = run_get_fisher_pvals(As,constants);
                end
                subspace(i).indices = get_max_edges(1-params.d_pval,subspace(i).size);
        end
    end
        
end