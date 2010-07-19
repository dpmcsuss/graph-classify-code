function plot_Lhats(Lhats,alg)

figure(4), clf, hold all

fs=12;      

colors{1}=0.75*[1 1 1];
colors{2}='k';

ls{1}='--';
ls{2}='-';
ls{3}=':';
ls{4}='-.';

fields=fieldnames(Lhats);

for i=1:alg.num_splits
    for j=1:alg.num_repeats
        for k=1:length(fields)
           temp.(fields{k})(i,j)=Lhats(i,j).(fields{k}); 
        end
    end
end

for i=1:length(fields)
    errorbar(alg.num_train_samples,mean(temp.(fields{i}),2),std(temp.(fields{i}),[],2),...
        'color',colors{mod(i,2)+1},'linestyle',ls{mod(i,3)+1},'linewidth',2);
end

axis([0 max(alg.num_train_samples)+1 0 0.5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)

legend(gca,fields,'Location','Best')
set(gca,'XTick',alg.num_train_samples,'XTickLabel',alg.num_train_samples)

if alg.save
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_Lhats'];
    print('-dpdf',figname)
    saveas(gcf,figname)
end
