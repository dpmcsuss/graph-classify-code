clear, clc
datadir = '/Users/jovo/Research/data/MRI/BLSA/will/2010_09_16/';
files   = dir(datadir); files(1:3)=[];
groups  = importdata('/Users/jovo/Research/data/MRI/BLSA/will/2010_04_14/BLSA32.txt');
s       = length(files)-3;
n       = 70;
As      = zeros(n,n,s);
ys      = zeros(1,s);
ids     = cell(1,s);
vm_slopes= zeros(1,s);

j=0;
for i=1:length(files)
    if strcmp(files(i).name(end-2:end),'csv')
        j=j+1;
        A= importdata([datadir files(i).name]);
        A(isnan(A))=0;
        A=A+A';
        As(:,:,j)=A(2:n+1,2:n+1);
        labels(j).id = files(i).name(16:17);

        for k=1:length(groups)
            if strcmpi(labels(j).id,groups{k}(1:2))
                if strcmp(groups{k}(12),'m')
                    labels(j).gender=1;
                else
                    labels(j).gender=0;
                end
                
                if strcmp(groups{k}(end),'t')
                    labels(j).hand=1;
                else
                    labels(j).hand=0;
                end
                
            end
        end

        %         for k=1:length(groups.textdata)
        %             if strcmpi(ids{j},groups.textdata{k,1})
        %                 m=k-1;
        %             end
        %         end
        %         ys(j) = groups.data(m,1);
        %         vm_slopes(j) = groups.data(m,2);
    end
end

% G = get_constants(As,ys);
% clear A datadir files groups i j k m n s
% save('~/Research/data/MRI/BLSA/will/BLSA42_data2','As','labels')

%% prepare data
% clear, clc, clf
% load('~/Research/data/MRI/BLSA/will/BLSA42')

for i=1:length(labels)
    A=As(:,:,i);
    Y = 0.2;                    % binarize 
    %Y = quantile(A(:),0.95);
    A(A<=Y)=0;
    As(:,:,i)=A;
    class_labels(i)=labels(i).gender;     % get Y from labels
end

adjacency_matrices=As;


alg.datadir = '~/Research/data/MRI/BLSA/will/';
alg.figdir  = '~/Research/figs/MRI/BLSA/';
alg.fname   = 'BLSA42';
alg.save    = 1;

alg.loopy   = 0;
alg.directed= 0;
alg.er      = 0;

% load([alg.datadir alg.fname])


% set algorithm parameters

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code

alg.ind_edge        = true;
alg.nb_ind          = find(triu(ones(constants.n),1));         % these graphs are undirected
alg.num_inc_edges   = 100; 
alg.num_coh_vertices= 10; 

alg.knn             = false;
alg.knn_vanilla     = true;

% test using in-sample training data
constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
[Lhatin Lvarin ind Pin yhatin] = graph_classify_ind_edge(adjacency_matrices,constants,alg); % compute in-sample classification accuracy
disp(Lhatin)

% test using hold-out training data

alg.num_splits=10;   % # of times to repeat
alg.num_repeats=10;  % # of times to repeat

%%
[Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg);

% make plots

constants   = get_constants(adjacency_matrices,class_labels);   % get constants like # edges, etc.
est_params  = get_params(adjacency_matrices,constants);         % estimate parameters from data

% plot_params(est_params,alg)                                     % plot params and estimated params
% plot_recovered_subspaces(constants,est_params,alg)              % plot recovered subspaces
[mean_hat std_hat] = plot_Lhats(Lhats,alg);                                           % plot misclassification rates
