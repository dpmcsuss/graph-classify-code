function xxx = plot_unlabeled_rates(Lhat_permuted,Lhats,Lhat_labeled,alg)

figure(5), clf, hold all

xxx.nb(1)=Lhat_permuted.nb;
xxx.tru(1)=Lhat_permuted.tru;
xxx.inc(1)=Lhat_permuted.inc;
xxx.coh(1)=Lhat_permuted.coh;

for j=1:alg.fw_max_iter
    xxx.nb(j+1)=Lhats{j}.nb;
    xxx.tru(j+1)=Lhats{j}.tru;
    xxx.inc(j+1)=Lhats{j}.inc;
    xxx.coh(j+1)=Lhats{j}.coh;
end

xxx.nb(j+2)=Lhat_labeled.nb;
xxx.tru(j+2)=Lhat_labeled.tru;
xxx.inc(j+2)=Lhat_labeled.inc;
xxx.coh(j+2)=Lhat_labeled.coh;


plot(0:alg.fw_max_iter+1,xxx.nb,'-','color',.75*[1 1 1],'linewidth',2)
plot(0:alg.fw_max_iter+1,xxx.tru,'-k','linewidth',2)
% plot(0:alg.fw_max_iter+1,xxx.inc,'-.k','linewidth',2)
% plot(0:alg.fw_max_iter+1,xxx.coh,'-.','color',0.75*[1 1 1],'linewidth',2)

axis('tight')
legend('nb','tru','inc','coh','location','NorthEast')
ylabel('misclassification rate')
xlabel('# of iterations')


if alg.save
    wh=[4 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_errorrates'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end