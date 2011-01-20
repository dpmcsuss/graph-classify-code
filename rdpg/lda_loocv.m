function [ Lhat ] = lda_loocv( features, classes, discrim, whiten )
%lda_loocv Does LOOCV and reports performance for discrim classifiers
%   Features is a dxn matrix, classes an 1xn vector of class labels and
%   discrim is a struct containing which classifiers to build: lda, dLda,
%   qda. Returns Lhat which is a nx1 struct array with a 1 for each
%   incorrectly classified sample. 
if nargin == 3
    whiten = false;
end

n = length(classes);
Lhat(n) = discrim;
Lsem = Lhat;
all_ind = struct('ytrn', 1:n, ...
                 'y0trn', find(classes == 0),...
                 'y1trn', find(classes == 1));
ind = all_ind;

            
for k=1:n
    ind.ytrn = all_ind.ytrn(all_ind.ytrn~=k);
    ind.y0trn = all_ind.y0trn(all_ind.y0trn~=k);
    ind.y1trn = all_ind.y1trn(all_ind.y1trn~=k);
    
    if ~whiten
        params = get_discriminant_params(features, ind, discrim);
        [Lhat(k) Lsem(k)] = discriminant_classifiers(features(:,k),classes(k),...
                                                params,discrim);
    else
%         [whiteFeat, wMat, uwMat] = fastica(features(:,ind.ytrn),...
%                                'only','white','verbose','off');
%         whiteFeat = wMat*features;
        
        whiteFeat = features-repmat(mean(features(:,ind.ytrn),2),1,n);
        whiteFeat = whiteFeat-repmat(mean(whiteFeat(:,ind.ytrn),2),1,n);
        [E,D] = eig(cov(whiteFeat(:,ind.ytrn)'),'nobalance');
        dim_keep = find(diag(D)>1e-7);
        E = E(:,dim_keep);
        %whiteFeat = whiteFeat(dim_keep,:);
        D = diag(diag(D(dim_keep,dim_keep)).^-.5);
        whiteFeat = E(dim_keep,:)*D*E'*whiteFeat;
        
%         end
        params = get_discriminant_params(whiteFeat, ind, discrim);
        [Lhat(k) Lsem(k)] = discriminant_classifiers(whiteFeat(:,k),...
                                              classes(k),params,discrim);
    end
                                         
end

end