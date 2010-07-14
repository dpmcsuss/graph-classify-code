clf, hold all
Ztrue=zeros(constants.n);
Ztrue(params.signal_subgraph_ind)=1;
for k=1:alg.num_splits
    jmax=length(inds{k});
    for j=1:jmax
        Zmat = zeros(constants.n);
        Zmat(inds{k}(j).inc) = 1;
        Ncorrect.inc(j) = sum(sum(Zmat(params.signal_subgraph_ind)));
        Nincorrect.inc(j) = sum(sum(Zmat~=Ztrue));
        
        Zmat = zeros(constants.n);
        Zmat(inds{k}(j).coh) = 1;
        Ncorrect.coh(j) = sum(sum(Zmat(params.signal_subgraph_ind)));
        Nincorrect.coh(j) = sum(sum(Zmat~=Ztrue));
    end
    true_positive_rate.inc(k) = mean(Ncorrect.inc)/params.num_signal_vertices^2;
    true_negative_rate.inc(k) = 1- mean(Nincorrect.inc)/((constants.n-params.num_signal_vertices)^2);
    std_correct.inc(k) = std(Ncorrect.inc)/params.num_signal_vertices^2;

    true_positive_rate.coh(k) = mean(Ncorrect.coh)/params.num_signal_vertices^2;
    true_negative_rate.coh(k) = 1- mean(Nincorrect.coh)/((constants.n-params.num_signal_vertices)^2);
    std_correct.coh(k) = std(Ncorrect.coh)/params.num_signal_vertices^2;
end

plot(true_negative_rate.inc,true_positive_rate.inc,'color',color1,'linestyle','+','linewidth',lw)
plot(true_negative_rate.coh,true_positive_rate.coh,'color',color2,'linestyle','o','linewidth',lw)

axis([0.97 1 0 1])
% set(gca,'XTick',1:k,'XTickLabel',2*alg.num_train_samples)
ylabel('true positive rate','fontsize',fs)
xlabel('true negative rate','fontsize',fs)
