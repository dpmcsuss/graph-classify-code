%% for cross-validation
alg.num_repeats = 1;        % # of times to repeat using the same # of hold-outs
alg.num_splits  = 1;        % # of different hold-outs
alg.naive_bayes = true;


%% for all sims
for i=1:length(sim)
    
    fprintf('\n\n**********************************')
    fprintf(['\nsimulation ' sim{i}.name '\n'])
    fprintf('**********************************\n\n')
    
    % use true parameters to get subspaces
    subspace.naive_bayes = find(sum(sim{i}.As,3)>0);
    subspace.incoherent  = sim{i}.model.m^2;
    subspace.coherent    = sim{i}.model.m;
    subspace = get_subspace(sim{i}.params,subspace);
    
    % use true parameters to estimate Lhat
    [Lhat_true Lsem_true] = naive_bayes_classify(sim{i}.As,sim{i}.targs,sim{i}.params,subspace);
    
    % check that using true parameters get Lstar, when Lstar is available
    if isfield(sim{i},'Lstar'),
        fprintf('L*      Lhat (sem)\n')
        fprintf('%1.2g \\in %1.2g (%1.2g) \n\n', sim{i}.Lstar,Lhat_true.naive_bayes,Lsem_true.naive_bayes)
        
        if Lhat_true.naive_bayes+2*Lsem_true.naive_bayes>=sim{i}.Lstar-1e-4
            disp('naive bayes classifier achieves Lhat + s.e. >= L*')
        else
            error('naive bayes classifier using true parameters fails to achieve Lhat within 1 s.e. of L*')
        end
        
        if Lhat_true.naive_bayes-2*Lsem_true.naive_bayes<=sim{i}.Lstar+1e-4
            disp('naive bayes classifier achieves Lhat - s.e. <= L*')
        else
            error('naive bayes classifier using true parameters fails to achieve Lhat within 1 s.e. of L*')
        end
        
    end
    
    % check that incoherent and coherent subspaces yield the same result
    % when using the true parameters and the subspace really is coherent
    if any(i==1:4)
        if sort(subspace.coherent)==sort(subspace.incoherent)
            disp('coherent and incoherent subspaces are the same, as they should be')
        else
            error('coherent and incoherent are different, but they should be the same')
        end
        
        if Lhat_true.coherent-Lhat_true.incoherent<1e-4
            disp('coherent and incoherent classifiers are the same, as they should be')
        else
            error('coherent and incoherent classifiers are different, but they should be the same')
        end
    end
    
    % get in-sample and hold-out performance
    constants = get_constants(sim{i}.As,sim{i}.targs);     % get constants to ease classification code
    
    % use true parameters to get subspaces
    subspace.naive_bayes = find(sum(sim{i}.As,3)>0);
    subspace.incoherent  = sim{i}.model.m^2;
    subspace.coherent    = sim{i}.model.m;
    fn = fieldnames(subspace);              % names of subspaces
    n_subspace = length(fn);                % # of different subspace
    
    % get insample performance
    [Lhat_in params_in subspace_in Lsem_in] = run_naive_bayes(sim{i}.As,constants,sim{i}.As,constants,subspace);
    
    % cross validate
    alg.s0_trn=sim{i}.model.s/4;     % list of #'s of hold-out samples for class 0
    alg.s1_trn=sim{i}.model.s/4;     % list of #'s of hold-out samples for class 1
    alg.num_train_samples=alg.s0_trn+alg.s1_trn;    % total # of training samples
    Lhat_holdout = classify_holdout(sim{i}.As,constants,alg,subspace);
    
    % check that in-sample outperforms holdout
    for j=1:n_subspace
        fprintf('\nLhat_in - 2*Lsem_in <= Lhat_holdout\n')
        fprintf('%1.2g (%1.2g)  < %1.2g\n', Lhat_in.(fn{j}), Lsem_in.(fn{j}), Lhat_holdout.(fn{j}))
        if Lhat_in.(fn{j})-2*Lsem_in.(fn{j})<=Lhat_holdout.(fn{j})
            disp([(fn{j}) ' in-sample performs better than hold-out'])
        else
            error([(fn{j}) ' hold-out outperforms in-sample'])
        end
    end
end
