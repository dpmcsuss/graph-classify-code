function plot_edge_identification_rates(inds,constants,alg)
%% plot misclassification
figure(5), clf, hold all

gray=0.75*[1 1 1];
lw=2;
fs=12;

for k=1:length(alg.num_train_samples)
    jmax=length(inds{k});
    for j=1:jmax
        Zmat = zeros(constants.n);
        Zmat(inds{k}(j).inc) = 1;
        Ncorrect.inc(j) = sum(sum(Zmat(alg.signal_subgraph_ind)));

        Zmat = zeros(constants.n);
        Zmat(inds{k}(j).coh) = 1;
        Ncorrect.coh(j) = sum(sum(Zmat(alg.signal_subgraph_ind)));
    end
    mean_correct.inc(k) = mean(Ncorrect.inc)/alg.num_signal_edges;
    std_correct.inc(k) = std(Ncorrect.inc)/alg.num_signal_edges;

    mean_correct.coh(k) = mean(Ncorrect.coh)/alg.num_signal_edges;
    std_correct.coh(k) = std(Ncorrect.coh)/alg.num_signal_edges;
end

errorbar(mean_correct.inc,std_correct.inc,'color','k','linestyle','-','linewidth',lw)
errorbar(mean_correct.coh,std_correct.coh,'color',gray,'linestyle','-','linewidth',lw)

axis([0 k 0 1])
set(gca,'XTick',1:k,'XTickLabel',alg.num_train_samples)
ylabel('edge detection rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)


if alg.save
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end