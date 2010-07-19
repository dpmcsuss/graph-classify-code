clear; clc

%% generate adjacency matrices and labels
female_left=xlsread('~/Research/necog/data/huiguang/Female_Left_ROIThick.xls');
female_right=xlsread('~/Research/necog/data/huiguang/Female_Right_ROIThick.xls');
female = [female_left female_right];
female(1,:)=[];
siz = size(female);
A0=zeros(siz(2),siz(2),siz(1));

for i=1:siz(1)
    A0(:,:,i)=female(i,:)'*female(i,:);
    for j=1:siz(2)
        for k=1:siz(2)
            A0(j,k,i)=A0(j,k,i)/(female(i,j)^2+female(i,k)^2);
        end
    end
end

male_left=xlsread('~/Research/necog/data/huiguang/male_Left_ROIThick.xls');
male_right=xlsread('~/Research/necog/data/huiguang/male_Right_ROIThick.xls');
male = [male_left male_right];
male(1,:)=[];

siz = size(male);
A1=zeros(siz(2),siz(2),siz(1));
for i=1:siz(1)
    A1(:,:,i)=male(i,:)'*male(i,:);
    for j=1:siz(2)
        for k=1:siz(2)
            A1(j,k,i)=A1(j,k,i)/(female(i,j)^2+female(i,k)^2);
        end
    end
end

siz0=size(A0);
siz1=size(A1);
adjacency_matrices(:,:,1:siz0(3))=A0;
adjacency_matrices(:,:,siz0(3)+1:siz0(3)+siz1(3))=A1;
class_labels=[zeros(1,siz0(3)) ones(1,siz1(3))];

for i=1:length(adjacency_matrices)
    A=adjacency_matrices(:,:,i);
    thr = quantile(A(:),0.5);
    A(A<=thr)=0;
    A(A>thr)=1;
    A=triu(A,+1);
    adjacency_matrices(:,:,i)=A;
end

alg.datadir = '~/Research/necog/data/huiguang/';
alg.figdir  = '~/Research/necog/figs/huiguang/';
alg.fname   = 'gender';
alg.save = 1;

if alg.save, save([alg.datadir alg.fname],'adjacency_matrices','class_labels'); end

%% set algorithm parameters

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

alg.ind_edge        = true;
alg.nb_ind          = find(triu(ones(constants.n),1));         % these graphs are undirected
alg.num_inc_edges   = 100; 
alg.num_coh_vertices= 10; 

alg.knn             = true;
alg.knn_vanilla     = true;

%% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

%% test using hold-out training data

alg.num_splits=10;   % # of times to repeat
alg.num_repeats=10;  % # of times to repeat

[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

%% make plots

constants   = get_constants(adjacency_matrices,class_labels);   % get constants like # edges, etc.
est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

plot_params(est_params,alg)                                     % plot params and estimated params
plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces
plot_Lhats(Lhats,alg)                                           % plot misclassification rates
