function [Lhats inds] = wrapper_hold_out_balanced_training_data(adjacency_matrices,class_labels,alg)

constants = get_constants(adjacency_matrices,class_labels);     % get constants to ease classification code
nmin=min(constants.s0,constants.s1);
inds{alg.num_splits} = [];

if isfield(alg,'signal_subgraph_ind'), Lhats.tru   = zeros(alg.num_splits,alg.num_repeats); end
Lhats.nb    = zeros(alg.num_splits,alg.num_repeats);
Lhats.inc   = zeros(alg.num_splits,alg.num_repeats);
Lhats.coh   = zeros(alg.num_splits,alg.num_repeats);

for i=1:alg.num_splits

    for j=1:alg.num_repeats

        ind0  = randperm(nmin);
        ind1  = randperm(nmin);

        y0trn = constants.y0(ind0(1:alg.num_train_samples(i)));
        y0tst = constants.y0(ind0(alg.num_train_samples(i)+1:alg.num_train_samples(i)+alg.num_test_samples));

        y1trn = constants.y1(ind1(1:alg.num_train_samples(i)));
        y1tst = constants.y1(ind1(alg.num_train_samples(i)+1:alg.num_train_samples(i)+alg.num_test_samples));

        Atrn = adjacency_matrices(:,:,[y0trn y1trn]);
        Atst = adjacency_matrices(:,:,[y0tst y1tst]);

        ytrn = [zeros(1,alg.num_train_samples(i)) ones(1,alg.num_train_samples(i))];
        ytst = [zeros(1,length(y1tst)) ones(1,length(y1tst))];

        Gtrn = get_constants(Atrn,ytrn);
        Gtst = get_constants(Atst,ytst);

        [Lhat Lvar ind] = graph_classify_ie(Atrn,Gtrn,alg,Atst,Gtst);

        if isfield(alg,'signal_subgraph_ind'),
            Lhats.tru(i,j)  = Lhat.tru;
            if j==1, inds{i} = ind; end % just store indices from first run, ignoring nb (as that is all), needed for plotting
        end
        Lhats.nb(i,j)   = Lhat.nb;
        Lhats.inc(i,j)  = Lhat.inc;
        Lhats.coh(i,j)  = Lhat.coh;

    end
end

if alg.save, save([alg.datadir alg.fname '_results'],'Lhats','inds'); end