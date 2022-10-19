function [PairTR] = Generate_Ant_index(numAnts)
k=1;
PairTR=zeros(numAnts*(numAnts-1)/2,2);
for i=1:numAnts
    for j=i+1:numAnts
        PairTR(k,1)=i;
    	PairTR(k,2)=j;
    	k=k+1;
    end
end
%load('././Forward_code/PairTR_0.mat');
%load('././Forward_code/PairTR.mat');
end