%% generate adjacency matrices and labels
clear; clc

load('~/Research/necog/data/huiguang/meanROIThick');
siz = size(meanROIThick);

nleft=20;
left=zeros(siz(2),siz(2),nleft);
for i=1:nleft
    left(:,:,i)=meanROIThick(i,:)'*meanROIThick(i,:);
    for j=1:siz(2)
        for k=1:siz(2)
            left(j,k,i)=left(j,k,i)/(meanROIThick(i,j)^2+meanROIThick(i,k)^2);
        end
    end
end

nright=21;
right=zeros(siz(2),siz(2),nright);
for i=1:nright
    right(:,:,i)=meanROIThick(i+nleft,:)'*meanROIThick(i+nleft,:);
    for j=1:siz(2)
        for k=1:siz(2)
            right(j,k,i)=right(j,k,i)/(meanROIThick(i+nleft,j)^2+meanROIThick(i+nleft,k)^2);
        end
    end
end

ncontrols=21;
control=zeros(siz(2),siz(2),ncontrols);
for i=1:ncontrols
    control(:,:,i)=meanROIThick(i+nleft+nright,:)'*meanROIThick(i+nleft+nright,:);
    for j=1:siz(2)
        for k=1:siz(2)
            control(j,k,i)=control(j,k,i)/(meanROIThick(i+nleft+nright,j)^2+meanROIThick(i+nleft+nright,k)^2);
        end
    end
end

ver='leftright_control';

switch ver
    case 'left_right'
        class0=left; class1=right;
    case 'left_control'
        class0=left; class1=control;
    case 'right_control'
        class0=right; class1=control;
    case 'leftright_control'
        class0(:,:,1:nleft)=left; class0(:,:,nleft+1:nleft+nright)=right; class1=control;
end
siz0=size(class0);
siz1=size(class1);

adjacency_matrices(:,:,1:siz0(3))=class0;
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=class1;
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];

for i=1:length(class_labels)
    A=adjacency_matrices(:,:,i);
    thr = quantile(A(:),0.5);
    A(A<=thr)=0;
    A(A>thr)=1;
    A=triu(A,+1);
    adjacency_matrices(:,:,i)=A;
end

alg.datadir = '~/Research/necog/data/huiguang/';
alg.figdir  = '~/Research/necog/figs/huiguang/';
alg.fname   = ['TLE_' ver];
alg.save = 1;

if alg.save, save([alg.datadir alg.fname],'adjacency_matrices','class_labels'); end

%% set algorithm parameters

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

alg.ind_edge        = true;
alg.nb_ind          = find(triu(ones(constants.n),1));         % these graphs are undirected
alg.num_inc_edges   = 100; 
alg.num_coh_vertices= 10; 

alg.knn             = false;
alg.knn_vanilla     = true;


%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% test using hold-out training data

alg.num_splits  = 5;    % # of folds in k-fold cross-validation
alg.num_repeats = 10;    % # of times to resample using each quadruple of s_trn0, s_trn1, s_tst0, s_tst1

[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

constants   = get_constants(adjacency_matrices,class_labels);   % get constants like # edges, etc.
est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

plot_params(est_params,alg)                                     % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces
plot_Lhats(Lhats,alg)                                           % plot misclassification rates
