function [ Lhat ] = lda_loocv( features, classes, discrim )
%lda_loocv Does LOOCV and reports performance for discrim classifiers
%   Features is a dxn matrix, classes an 1xn vector of class labels and
%   discrim is a struct containing which classifiers to build: lda, dLda,
%   qda. Returns Lhat which is a nx1 struct array with a 1 for each
%   incorrectly classified sample. 
n = length(classes);
Lhat(n) = discrim;
all_ind = struct('ytrn', 1:n, ...
                 'y0trn', find(classes == 0),...
                 'y1trn', find(classes == 1));
ind = all_ind;
            
for k=1:n
    ind.ytrn = all_ind.ytrn(all_ind.ytrn~=k);
    ind.y0trn = all_ind.y0trn(all_ind.y0trn~=k);
    ind.y1trn = all_ind.y1trn(all_ind.y1trn~=k);
    
    params = get_discriminant_params(features, ind, discrim);
    [Lhat(k) Lsem(k)] = discriminant_classifiers(features(:,k),classes(k),...
                                                params,discrim);
end

end