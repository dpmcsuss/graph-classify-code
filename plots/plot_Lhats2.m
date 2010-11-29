function [Lhat_avg Lhat_std] = plot_Lhats2(Lhats,alg)

figure(4), clf, hold all

fs=12;      

colors{1}=0.75*[1 1 1];
colors{2}='k';

ls{1}='--';
ls{2}='-';
ls{3}=':';
ls{4}='-.';

fn=fieldnames(Lhats);
s_trn=alg.num_train_samples;    % total # of training samples

for i=1:length(fn)
    Lhat_avg.(fn{i})=mean(Lhats.(fn{i}));
    Lhat_std.(fn{i})=std(Lhats.(fn{i}));
    
    errorbar(s_trn,Lhat_avg.(fn{i}),Lhat_std.(fn{i}),...
        'color',colors{mod(i,2)+1},'linestyle',ls{mod(i,3)+1},'linewidth',2);
end

axis([0 max(s_trn)+1 0 0.5])
ylabel('misclassification rate','fontsize',fs)
xlabel('# training samples','fontsize',fs)

for k=1:length(fn)
    underscore = strfind(fn{k},'_');
    if ~isempty(underscore)
        fn{k}(underscore)=' ';
    end
end

legend(gca,fn,'Location','SouthWest')
set(gca,'XTick',s_trn,'XTickLabel',s_trn)

if alg.save
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_Lhats'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
