function [Lhat_avg Lhat_sem] = plot_Lhats2(Lhats,alg)


fn=fieldnames(Lhats);

if isfield(alg,'num_train_samples'),
    s_trn=alg.num_train_samples;
elseif isfield(alg,'trn_s')
    s_trn=alg.trn_s;
else
    s_trn=1;
    display('std reported instead of sem because @ of training samples not specified')
end

plot_do=1; if isfield(alg,'plot'), if alg.plot==0; plot_do=0; end, end
save_do=0; if isfield(alg,'save'), if alg.save==0; save_do=0; end, end

if save_do && plot_do
    figure(4), clf, hold all
    fs=12;
    
    colors{1}=0.75*[1 1 1];
    colors{2}='k';
    
    ls{1}='--';
    ls{2}='-';
    ls{3}=':';
    ls{4}='-.';
end

for i=1:length(fn)
    Lhat_avg.(fn{i})=mean(Lhats.(fn{i}));
    Lhat_sem.(fn{i})=std(Lhats.(fn{i}))./sqrt(s_trn);
    
    if plot_do
        errorbar([0 s_trn],[0.5 Lhat_avg.(fn{i})],[0 Lhat_sem.(fn{i})],...
            'color',colors{mod(i,2)+1},'linestyle',ls{mod(i,3)+1},'linewidth',2);
    end
end


for k=1:length(fn)
    underscore = strfind(fn{k},'_');
    if ~isempty(underscore)
        fn{k}(underscore)=' ';
    end
end

if plot_do
    axis([0 max(s_trn)+1 0 0.5])
    ylabel('misclassification rate','fontsize',fs)
    xlabel('# training samples','fontsize',fs)
    
    legend(gca,fn,'Location','SouthWest')
    set(gca,'XTick',s_trn,'XTickLabel',s_trn)
end

if save_do && plot_do
    wh=[3 2]*1.5;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[alg.figdir alg.fname '_Lhats'];
    print('-dpdf',figname)
    print('-deps',figname)
    saveas(gcf,figname)
end
