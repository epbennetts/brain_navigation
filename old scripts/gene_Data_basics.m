%we have a matrix that's like 20,000x213
clear vars;
close all;
load("C:\Users\elobe\Documents\1_HONOURS\Gene\GeneDataFile.mat", '-mat')

% %range == between -1 and 1 
% in_range = 0;
% low_range = 0; 
% high_range = 0;
% other = 0;
% smallest = 0; 
% largest = 0; 

gene_expr = zeros(1,19419);
%loop over areas
for i = 1:5
    %loop over genes
    for j = 1:19419
        % how different is this gene in area X vs ~X? 
        %make array of first row
            if ~isnan(geneData(i,j)) 
            gene_expr(j) = geneData(i,j);
            elseif j ~= 1 %and obv isnan
                gene_expr(j) = gene_expr(j-1);
            %else it is NaN but is alsothe first element in the row 
            %in which case we just leave it at 0    
            end
        
%         %testing max and min values
%         if ZgeneData(i, j) < -1
%             low_range = low_range + 1;
%             if ZgeneData(i, j) < smallest
%                 smallest = ZgeneData(i, j);
%             end
%         elseif ZgeneData(i, j) > 1
%             high_range = high_range + 1;
%             if ZgeneData(i, j) > largest
%                 largest = ZgeneData(i, j);
%             end
%         elseif ZgeneData(i, j) <= 1 && ZgeneData(i, j) >= -1 
%             in_range = in_range + 1;  
%         else 
%             other = other + 1;
%         end
 
    end
    figure()
    plot(1:19419,gene_expr)
    xlabel('genes');
    ylabel('expression (z-scored twice?)');
    
end
%wow, I didn't even need a for loop...

