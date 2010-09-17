clear, clc
datadir = '/Users/joshyv/Research/data/MRI/BLSA/will/2010_08_16/';
files   = dir(datadir); files(1:3)=[];
groups  = importdata('/Users/joshyv/Research/data/MRI/BLSA/will/2010_04_14/BLSA32.txt');
s       = length(files)-3;
n       = 70;
As      = zeros(n,n,s);
ys      = zeros(1,s);
ids     = cell(1,s);
vm_slopes= zeros(1,s);

j=0;
for i=1:length(files)
    if strcmp(files(i).name(end-2:end),'csv')
        j=j+1;
        A= importdata([datadir files(i).name]);
        A(isnan(A))=0;
        A=A+A';
        As(:,:,j)=A(2:n+1,2:n+1);
        labels(j).id = files(i).name(16:17);

        for k=1:length(groups)
            if strcmpi(labels(j).id,groups{k}(1:2))
                if strcmp(groups{k}(12),'m')
                    labels(j).gender=1;
                else
                    labels(j).gender=0;
                end
                
                if strcmp(groups{k}(end),'t')
                    labels(j).hand=1;
                else
                    labels(j).hand=0;
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

% G = get_constants(As,ys);
% clear A datadir files groups i j k m n s
save('~/Research/data/MRI/BLSA/will/BLSA42','As','labels')
