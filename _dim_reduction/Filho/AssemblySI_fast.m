function [CorrMatrix,threshold]=AssemblySI_fast(Patterns1,Patterns2,n_aleat,p_lim)
% Code to calculate the cosine similarity between the patterns of 2
% assembly matrices and define its threshold.
% Inputs:
%   1) Patterns1: N x M matrix in which N comprise the number of neurons in
%   the dataset, and M, the number of assemblies.
%   2) Patterns2: N x P matrix in which N comprise the same number of 
%   neurons as Patterns 1, and P, the number of assemblies.
%   3) n_aleat: positive scalar with the number of aleatorizations to determine
%   significance (Default: 1000);
%   4) p_lim: scalar from 0 to 1 determining the p value considered as
%   limit for significance(Default: 0.001).
%
% Outputs:
%   1) SimMatrix: M x P matrix corresponding to the similarity index
%   calculated between the M assemblies of Patterns1 and the P assemblies
%   of Patterns2.
%   2) threshold: empirically calculated similarity index threshold for 
%   significance between assemblies. 
% Daniel Almeida Filho, almeidafilhodg@ucla.edu

if nargin<3
    n_aleat=1000;
    p_lim=0.001;
end


CorrMatrix=abs(corr(Patterns1,Patterns2));

NumberNeurons=size(Patterns1,1);

buffer=[];
fake1=[];
fake2=[];

for i=1:n_aleat
    fake1(:,:)=Patterns1(randperm(NumberNeurons),:);
    fake2(:,:)=Patterns2(randperm(NumberNeurons),:);
    temp = abs(corr(fake1,fake2));
    buffer = [buffer reshape(temp,1,[])];
end

r =sort(buffer,'descend');
threshold=r(round(p_lim*length(r)));
