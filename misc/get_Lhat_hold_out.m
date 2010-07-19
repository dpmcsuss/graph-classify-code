function [Lhats alg inds] = get_Lhat_hold_out(adjacency_matrices,class_labels,alg)
% this function computes trains and classifies using held-out data
% unbalanced training data is permitted
% number test samples can be a function of number of training samples
% 
% INPUT:
%   adjacency_matrices: 3rd-order tensor of adjacency_matrices
%   class_labels: vector of class labels
%   alg: specifies which algorithms to use, and which parameters for each alg
% 
% OUTPUT:
%   Lhats: structure with Lhat from each algorithm run
%   alg: some alg fields can be set by default in here, and must be output for displaying results
%   inds: indices of signal subgraphs found by various algorithms

constants = get_constants(adjacency_matrices,class_labels);         % get constants to ease classification code

if ~isfield(alg,'num_splits'), alg.num_splits = 1; end              % # of splits in k-fold cross-validation
if ~isfield(alg,'num_repeats'), alg.num_repeats = 1; end            % # of times to repeat each split

if ~isfield(alg,'num_class0_train_samples')                         % # of samples to train class 0 parameters per fold
    alg.num_class0_train_samples    = round(linspace(1,constants.s0-1,alg.num_splits));
end

if ~isfield(alg,'num_class0_test_samples')                          % # of samples to test class 1 per fold
    alg.num_class0_test_samples     = constants.s0-alg.num_class0_train_samples;
end

if ~isfield(alg,'num_class1_train_samples')                         % # of samples to train class 1 parameters per fold
    alg.num_class1_train_samples    = round(linspace(1,constants.s1-1,alg.num_splits));
end

if ~isfield(alg,'num_class1_test_samples')                          % # of samples to test class 1 per fold
    alg.num_class1_test_samples     = constants.s1-alg.num_class1_train_samples;
end

alg.num_train_samples   = alg.num_class0_train_samples+alg.num_class1_train_samples;    % total number of training samples per fold
alg.num_test_samples    = alg.num_class0_test_samples+alg.num_class1_test_samples;      % total number of testing samples per fold

if any(alg.num_class0_train_samples+alg.num_class0_test_samples>constants.s0) || ...    % check make sure there are not more samples to use in testing and training then we generated
        any(alg.num_class1_train_samples+alg.num_class1_test_samples>constants.s1),
    error('cannot have more testing and training samples than total samples');
end

if ~isfield(alg,'ind_edge'),    alg.ind_edge    = true;     end
if ~isfield(alg,'knn'),         alg.knn         = false;    end

inds{alg.num_splits} = [];                                       % initialize memory for indices (not really though)

if alg.ind_edge
    if isfield(alg,'signal_subgraph_ind'), Lhats.tru = zeros(alg.num_splits,alg.num_repeats); end
    Lhats.nb    = zeros(alg.num_splits,alg.num_repeats);
    Lhats.inc   = zeros(alg.num_splits,alg.num_repeats);
    Lhats.coh   = zeros(alg.num_splits,alg.num_repeats);
end

for i=1:alg.num_splits

    for j=1:alg.num_repeats

        ind0  = randperm(constants.s0);
        ind1  = randperm(constants.s1);

        y0trn = constants.y0(ind0(1:alg.num_class0_train_samples(i)));
        y0tst = constants.y0(ind0(alg.num_class0_train_samples(i)+1:alg.num_class0_train_samples(i)+alg.num_class0_test_samples(i)));

        y1trn = constants.y1(ind1(1:alg.num_class1_train_samples(i)));
        y1tst = constants.y1(ind1(alg.num_class1_train_samples(i)+1:alg.num_class1_train_samples(i)+alg.num_class1_test_samples(i)));

        Atrn = adjacency_matrices(:,:,[y0trn y1trn]);
        Atst = adjacency_matrices(:,:,[y0tst y1tst]);

        ytrn = [zeros(1,length(y0trn)) ones(1,length(y1trn))];
        ytst = [zeros(1,length(y0tst)) ones(1,length(y1tst))];

        Gtrn = get_constants(Atrn,ytrn);
        Gtst = get_constants(Atst,ytst);

        if alg.ind_edge
            [Lhat Lvar ind] = graph_classify_ind_edge(Atrn,Gtrn,alg,Atst,Gtst);
            if isfield(alg,'signal_subgraph_ind'),
                Lhats.tru(i,j)  = Lhat.tru;
                if j==1, inds{i} = ind; end % just store indices from first run, ignoring nb (as that is all), needed for plotting
            end
            Lhats.nb(i,j)   = Lhat.nb;
            Lhats.inc(i,j)  = Lhat.inc;
            Lhats.coh(i,j)  = Lhat.coh;
        end

        if alg.knn
            Lhat_knn = graph_classify_knn(Atrn,Gtrn,alg,Atst,Gtst);
            Lhats.knn(i,j)=Lhat_knn.knn;
        end

    end
end

if alg.save, save([alg.datadir alg.fname '_results'],'Lhats','inds'); end