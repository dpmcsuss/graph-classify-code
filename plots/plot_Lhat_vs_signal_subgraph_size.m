colors{3}='k';
colors{1}=0.75*[1 1 1];
colors{4}='k';
colors{2}=0.75*[1 1 1];

ls{1}='--';
ls{2}='-';
ls{3}='-';
ls{4}='--';

for i=1:70, blowup(i)=choose(i,2); end

figure(2), clf, hold on
for i=1:n_subspaces
    if ~strcmp(fn{i},'coherent')
        plot(num_inc_edges,Lhats.(fn{i}),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{i},'linewidth',2);
    else
        plot(blowup(7:end),Lhats.(fn{i})(7:end),... Lsems.(fn{i}),...
            'color',colors{i},'linestyle',ls{i},'linewidth',2);
    end
end

axis('tight')
legend('incoherent','inc., scaled','fisher','coherent')
set(gca,'XScale','log')
set(gca,'XLim',[20 2415])
set(gca,'YLim',[.37 .5])

xlabel('size of signal subgraph','fontsize',12)
ylabel('misclassification rate','fontsize',12)

if alg.save
    wh=[6 3];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_Lhat_vs_m'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
