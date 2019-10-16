function [gene_select_result, gene_number, NSLRG_descending_order,select_number,der2_value]=NSLRG_S(X,Z_hat)
%% %% this matlab code is used to NSLRG and Score function and obtain  the selected feature genes.
%% %%----------------------------
%% %%----------------------------
% inputs
%       X -- data matrix, M*N (M>N)
%       Z_hat -- the low-rank symmetric matrix., N*N

% output
%        gene_select_result : Selected subset of genes
%        gene_number: The line number of the selected gene in the original data
%        NSLRG_descending_order: The scores are sorted in descending order
%        select_number: The number of genes selected
%        der2_value :Fractional curve fitting value

%%  %%

%% %% main code

[m,n] = size(X);

W=Z_hat;
NSLRG=[];% score
gene_sample_similarity=[]; %The sum of the similarities of a sample to all other samples in one gene.
sum_gene_sample_similarity=[]; %The sum of the similarities of all samples to one another in a given gene.
gene_VAR=(std(X,0,2)).^2; %Find the variance of expression data of a gene.

for r=1:1:m  
    for i=1:1:n
        for j=1:1:n
        gene_sample_similarity(j,1)=((X(r,i)-X(r,j))^2)*W(i,j);
        end
        sum_gene_sample_similarity(i)=sum(gene_sample_similarity);
    end
    NSLRG(r,1)=sum(sum_gene_sample_similarity)/gene_VAR(r);
end


NSLRG_descending_order=sort(NSLRG(:)); %Ascending order points

%Determination of turning point of curve
data=NSLRG_descending_order';
y=data(isfinite(data));
y=y(isfinite(y));
x=1:length(y); 

%fitting
p=polyfit(x,y,30); 
%der
yy=polyder(p);
yyy=polyder(yy);
der2_value=polyval(yyy,x);

figure
subplot(1,2,1);
plot(x,y);
hold on
subplot(1,2,1);
plot(x,polyval(p,x));
hold on
subplot(1,2,2);
plot(x,der2_value);
drawnow;
select_number=input('input the number of select feature base on the figure:');
 
close;



NSLRG_select=NSLRG_descending_order(1:1:select_number,1); % The first select_number score.
gene_number=zeros(select_number,1); % A subset of the original data, X, put the selected gene here.

FIND=unique(NSLRG_select);
[select_number_unique ,~]=size(FIND);
number=1;
for number_1=1:1:select_number_unique %Select the line number of the original genetic data corresponding to a score.
    gene_find_number=find(NSLRG==FIND(number_1)); 
    
    % Determine if there are genes with the same score, and if there are genes with the same score, select their line Numbers
    [number_m,~] = size(gene_find_number);  
    if number_m==1
    gene_number(number)=gene_find_number;
     number=number+1;
    else
        for loop=1:1:number_m
        gene_number(number)=gene_find_number(loop,1);
        number=number+1;
        end
    end
end

%The genetic data of the corresponding line number were extracted
gene_select_result=X(gene_number,:);


end









