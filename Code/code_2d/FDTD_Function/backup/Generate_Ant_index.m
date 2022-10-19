function [PairTR] = Generate_Ant_index(numAnts)
k=1;
for i=1:numAnts
    for j=1:numAnts
        if i<j
            PairTR(k,1)=i;
            PairTR(k,2)=j;
            k=k+1;
        end
    end
end
end