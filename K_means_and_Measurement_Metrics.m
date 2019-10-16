  function [total_label,total_res,result ] = K_means_and_Measurement_Metrics (Samplecategory,kmeans_number,kmeans_data)
%% %%
%   input:
%       Samplecategory: label
%       kmeans_number: run number
%       gene_select_result
%   output:
%       total_label kmeans Result label
%       total_res   the mapping after the corresponding tag
%       result.acc:  [maximum, minimum, average, variance,  average,  variance]
%       result.ri_1 Same as above
%       result.REC  Same as above
%       result.PRE  Same as above
%       result.F1  Same as above
 
%% %% 
%Kmeans
gnd=Samplecategory;
k=numel(unique(gnd));
total_label=[];
total_res=[];
ACC=[];
ri_1=[];
ri_2=[];
NMI=[];
recall=[];
precision=[];
F1_Measure=[];

for clustering_number=1:1:kmeans_number
    label=kmeans(kmeans_data,k);
    total_label(end+1,:)=label;
%% %% ACC
    res=bestMap(gnd,label);
    total_res(end+1,:)=res;
    ACC(end+1)=length(find(gnd==res))/length(gnd);
%% %%rand index
    p1=gnd;
    p2=label;
    ri_1(end+1)=rand_index(p1, p2);
    ri_2(end+1)=rand_index(p1, p2, 'adjusted');
 
%% %% [REC,PRE,F1]
    [REC,PRE,F1]=REC_PRE_F1(Samplecategory,res);
    recall(end+1,:)=REC;
    precision(end+1,:)=PRE;
    F1_Measure(end+1,:)=F1;
end
%% Accuracy of the results
MAX_ACC=max(ACC);% Maximum
MIN_ACC=min(ACC);% minimum
MEAN_ACC=sum(ACC)/clustering_number; %average accuracy
VAR_ACC=var(ACC);%The variance

varcut_ac=[];
cut_acc=0;
cut_acc_n=0;
for i=1:1:clustering_number
    if ACC(i)>0.7
        cut_acc=cut_acc+ACC(i);
        cut_acc_n=cut_acc_n+1;
        varcut_ac(cut_acc_n)=ACC(i);
    end
end
cut_ac=cut_acc/cut_acc_n;
VAR_cut_ac=var(varcut_ac);


result.acc=[MAX_ACC,MIN_ACC,MEAN_ACC,VAR_ACC, cut_ac,VAR_cut_ac];
% [maximum accuracy, minimum accuracy, average accuracy, variance]

%% rand index ri_1
MAX_ri_1=max(ri_1);% Maximum
MIN_ri_1=min(ri_1);% minimum
MEAN_ri_1=sum(ri_1)/clustering_number; %average
VAR_ri_1=var(ri_1); %The variance

varcut_ri_11=[];
cut_ri_1=0;
cut_ri_1_n=0;
for i=1:1:clustering_number
    if ri_1(i)>0.7
        cut_ri_1=cut_ri_1+ri_1(i);
        cut_ri_1_n=cut_ri_1_n+1;
         varcut_ri_11( cut_ri_1_n)=ri_1(i);
    end
end
cut_ri_11=cut_ri_1/cut_ri_1_n;
VAR_ri_11=var(varcut_ri_11);

result.ri_1=[MAX_ri_1,MIN_ri_1,MEAN_ri_1,VAR_ri_1, cut_ri_11,VAR_ri_11];
% [maximum rand index, minimum rand index, average rand index, variance]

%% %% [REC,PRE,F1]
%recall
MAX_REC=max(recall);% Maximum
MIN_REC=min(recall);% minimum
MEAN_REC=sum(recall)/clustering_number; %average
VAR_REC=var(recall);
result.REC=[MAX_REC,MIN_REC,MEAN_REC,VAR_REC];
% [maximum recall, minimum recall, average recall, variance]

%precision
MAX_PRE=max(precision);% Maximum
MIN_PRE=min(precision);% minimum
MEAN_PRE=sum(precision)/clustering_number; %average
VAR_PRE=var(precision);
result.PRE=[MAX_PRE,MIN_PRE,MEAN_PRE,VAR_PRE];
%  [maximum precision, minimum precision, average precision, variance]

%F1_Measure
MAX_F1 = max(F1_Measure);% Maximum
MIN_F1=min(F1_Measure);% minimum
MEAN_F1=sum(F1_Measure)/clustering_number; %average
VAR_F1=var(F1_Measure);
result.F1=[MAX_F1,MIN_F1,MEAN_F1,VAR_F1]; 
%  [maximum F1, minimum , average , variance]


