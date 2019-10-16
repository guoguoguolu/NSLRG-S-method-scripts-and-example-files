function [REC,PRE,F1]=REC_PRE_F1(Samplecategory,res)
%% %%%%
%   input  Samplecategory ture label
%           label    predivting laber
%   output  recall   
%               precision
%               F1_Measure
%% %%
 

%% %% main code
m=length(unique(Samplecategory));
    REC=[];
    PRE=[];
    F1=[];
    for i=1:1:m
        
    Orgi=length(find(Samplecategory==i));
    Pred=length(find(res==i));
    
    T_Pred=length(intersect(find(Samplecategory==i),find(res==i)));
    
    precision=T_Pred/Pred;
    recall=T_Pred/Orgi;
    
    F1_Measure=2/(1/precision + 1/recall);  
    
    REC(1,i)=recall;
    PRE(1,i)=precision;
    F1(1,i)=F1_Measure;
    end








