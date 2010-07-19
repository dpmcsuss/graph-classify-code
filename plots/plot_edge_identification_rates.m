function plot_edge_identification_rates(inds,constants,alg)
%% plot misclassification
figure(5), clf, hold all

gray=0.75*[1 1 1];
lw=2;
fs=12;

for i=1:alg.num_splits
    for j=1:alg.num_repeats
        kmax=length(inds{i,j});
        for k=1:kmax
            Zmat = zeros(constants.n);
            Zmat(inds{i,j}(k).inc) = 1;
            Ncorrect.inc(k) = sum(sum(Zmat(alg.signal_subgraph_ind)));

            Zmat = zeros(constants.n);
            Zmat(inds{i,j}(k).coh) = 1;
            Ncorrect.coh(k) = sum(sum(Zmat(alg.signal_subgraph_ind)));
        end
        mean_correct.inc(j) = mean(Ncorrect.inc)/alg.num_signal_edges;
        mean_correct.coh(j) = mean(Ncorrect.coh)/alg.num_signal_edges;
    end
    avg_correct.inc(i)=mean(mean_correct.inc);
    avg_correct.coh(i)=mean(mean_correct.coh);

    ste_correct.inc(i)=std(mean_correct.inc);
    ste_correct.coh(i)=std(mean_correct.coh);
    
end
errorbar(avg_correct.inc,ste_correct.inc,'color','k','linestyle','-','linewidth',lw)
errorbar(avg_correct.coh,ste_correct.coh,'color',gray,'linestyle','-','linewidth',lw)

axis([0 alg.num_splits 0 1])
set(gca,'XTick',1:k,'XTickLabel',alg.num_train_samples)
ylabel('edge detection rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)
legend('inc','coh')

if alg.save
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end