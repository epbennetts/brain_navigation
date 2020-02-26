%we have a matrix that's like 20,000x213
clear vars;
load("C:\Users\elobe\Documents\1_HONOURS\Gene\GeneDataFile.mat", '-mat')

ZgeneData = zscore(geneData);

% %range == between -1 and 1 
% in_range = 0;
% low_range = 0; 
% high_range = 0;
% other = 0;
% 
% smallest = 0; 
% largest = 0; 

gene_expr = zeros(1,19419);
j = 1:19419;
%loop over areas
for i = 1
    figure(i)
    plot(j,geneData(i,:)) 
end