clear vars;
close all;
load("C:\Users\elobe\Documents\1_HONOURS\Gene\GeneDataFile.mat", '-mat')

ZgeneData = zscore(geneData);

row_entropy = zeros(213);

%loop over areas
for i = 1:5
    gene_expr = zeros(1,19419);
    gene_entropy = zeros(1,19419);
    
    %entropy of dataset
    S = -(1/213*log2(1/213) + 212/213*log2(212/213));
    
    %loop over genes
    for j = 1:19419
        % how different is this gene expressed in area X vs ~X? 
        %make arrays of first rows
            if isnan(ZgeneData(i,j)) 
                ZgeneData(i,j) = 0;
            end
            gene_expr(j) = ZgeneData(i,j);
           % gene_entropy(j) = 
            %else it is NaN
            %in which case we just leave it at 0    
            
    end
    figure()
    plot(1:19419,gene_expr)
    %should be plotting expression of one gene in all areas 
    %instead of all genes in 1 area

end
