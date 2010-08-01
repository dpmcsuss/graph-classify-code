clear, clc
datadir = '~/Research/necog/data/BLSA/will/BLSA82/';
files   = dir(datadir);
groups  = importdata('~/Research/necog/data/BLSA/will/BLSA82/BLSA32.txt');
s       = length(files)-3;
n       = 70;
As      = zeros(n,n,s);
ys      = zeros(1,s);
ids     = cell(1,s);
vm_slopes= zeros(1,s);

for i=3:length(files)
    j=i-2;
    if strcmp(files(i).name(end-2:end),'csv')
        A= importdata([datadir files(i).name]);
        A(isnan(A))=0;
        A=A+A';
        As(:,:,j)=A(2:n+1,2:n+1);
        ids{j} = files(i).name(16:17);

        for k=1:length(groups)
            if strcmpi(ids{j},groups{k}(1:2))
                if strcmp(groups{k}(12),'m')
                    gender(j)=1;
                else
                    gender(j)=0;
                end
                
                if strcmp(groups{k}(end),'t')
                    hand(j)=1;
                else
                    hand(j)=0;
                end
                
            end
        end

        %         for k=1:length(groups.textdata)
        %             if strcmpi(ids{j},groups.textdata{k,1})
        %                 m=k-1;
        %             end
        %         end
        %         ys(j) = groups.data(m,1);
        %         vm_slopes(j) = groups.data(m,2);
    end
end

G = get_constants(As,ys);
% clear A datadir files groups i j k m n s
save('~/Research/necog/data/BLSA/will')
