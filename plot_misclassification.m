function plot_misclassification(Lhats,inds,constants,alg)
%% plot misclassification
figure(2), clf, hold all

gray=0.75*[1 1 1];
lw=2;
fs=12;


if isfield(alg,'signal_subgraph_ind'), subplot(121), hold all; end
h(1)=errorbar(alg.num_train_samples,mean(Lhats.nb,2),var(Lhats.nb,[],2),'color',gray,'linestyle','--','linewidth',lw);
h(2)=errorbar(alg.num_train_samples+0.2,mean(Lhats.inc,2),var(Lhats.inc,[],2),'color','k','linestyle','-','linewidth',lw);
h(3)=errorbar(alg.num_train_samples+0.4,mean(Lhats.coh,2),var(Lhats.coh,[],2),'color',gray,'linestyle','-','linewidth',lw);
if isfield(alg,'signal_subgraph_ind'), h(4)=errorbar(alg.num_train_samples+0.6,mean(Lhats.tru,2),var(Lhats.tru,[],2),'color','k','linestyle','--','linewidth',lw); end

axis([0 max(alg.num_train_samples)+1 0 0.5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)
legend(gca,'naive Bayes','incoherent','coherent','optimal','Location','Best')
set(gca,'XTick',alg.num_train_samples,'XTickLabel',alg.num_train_samples)

if isfield(alg,'signal_subgraph_ind'), subplot(122), hold all;

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

end


if alg.save
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_rates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end
