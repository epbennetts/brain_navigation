    %TO DO
    %later: have clean genes matrix saved somewhere (instead of processing every time)
    %generalise to a function that uses n num of genes and just goes through
    %everything without plotting
    %generalise further to any area

%     numgenes = 1;
%     samples = 19114;
%     area = "Isocortex";
%     prev_best_genes = [1237, 7, 48];
    

function [indexOrder, accuracies_ranked, genenames_ranked, thresholds_all] = DT_classification_multiple(numgenes, samples, area, prev_best_genes)
    %Input: 
    %numgenes: Number of genes to use in each DT
    %samples: number of "samples" to complete the algorithm on (can choose samples < cols for
    %testing)
    %prev_best_genes: array of the indices of the best genes for each
    %iteration of increasing numgenes. i.e. for numgenes = 1, don't need
    %anything, for numgenes = 2, need indexOrder(1) of the previous iteration, 
    %for numgenes = 3, need array of size 2 with the 2 prev ones, etc.
    %ideally this would eventually be done recursively within the function
    %rather than manually outside
    %Function:
    %The function will train a tree with those genes and test it
    %against the real classes. Will also plot the accuracies histogram
    %Returns: 
    %indexOrder: gene index sorted by decreasing accuracy
    %accuracies_ranked: Accuracies from highest to lowest
    %genenames_ranked: Gene names ordered by highest to lowest accuracy
    %thresholds_all: DT thresholds for every tree


    [genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans(area);
    [rows, cols] = size(genes);

    %num of trees we want to plot (not in order)
    %plotting = 1:10;
    %plotting = [];

    %set vars and empty arrays
    accuracies_new = zeros(samples,1);
    conf_matrices = zeros(samples,4);
    thresholds_all = NaN(samples, numgenes);
    %genes_used = cell(samples, 2);

    try
        if (numgenes == 2)
            a = prev_best_genes(1);
        elseif (numgenes == 3)
            a = prev_best_genes(1);
            b = prev_best_genes(2);
        end
    catch
        disp("Not enough elements in prev_best_genes array");
    end


    %loop, make trees, classify and evaluate
    for i = 1:samples

        %set gene pair
        if numgenes == 1
            gene_combo = genes(:,i);
        elseif numgenes == 2
            gene_combo = [genes(:,a) genes(:,i)]; %gene = genes(:,i:i+1);
        elseif numgenes == 3
            gene_combo = [genes(:,a) genes(:,b) genes(:,i)];
        end

        %train
        tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);

        %check first predictor (should be best gene: x1)
        %doesn't exactly work as doesn't always show second split?
        %gene_used = tree.CutPredictor;
        %genes_used{i,1} = gene_used{1};
        %genes_used{i,2} = gene_used{2};

        %thresholds
        thresholds_raw = tree.CutPoint; %threshold with NaNs
        thresholds = [];
        %delete nans
        for j = 1:size(thresholds_raw,1)
            thr = thresholds_raw(j);
            if ~isnan(thr)
                thresholds = [thresholds thr];
                if j>1
                end
            end
        end
        %1st threshold
        if ~isnan(thresholds)
            t1 =  thresholds(1,1);
            thresholds_all(i,1) = t1;
        end
        %2nd & 3rd threshold
        if(size(thresholds,2) == 2)
            t2 = thresholds(1,2);
            thresholds_all(i,2) = t2;
        elseif(size(thresholds,2) == 3)
            t3 = thresholds(1,3);
            thresholds_all(i,3) = t3;
        else
            sprintf("There were more than 3 thresholds");
        end
        %gene_index_used = tree.CutPredictorIndex;

        %test
        labels = predict(tree, gene_combo);
        %threshold_old = edge(tree, gene_combo, classes);

       % calc old accuracy (for reference)
%        correct = zeros(rows,1);
%         for j = 1:rows
%             correct(j) = (strcmp(labels{j}, classes(j)));
%         end
%         accuracies_old(i) = mean(correct,1);

        %calc confusion matrix
        %targets
        labels_t = labels(targetIndices == 1);
        tp = size(labels_t(strcmp(labels_t,"target")),1);
        fn = size(labels_t(strcmp(labels_t,"~target")),1);
        %non-targets
        labels_nt = labels(targetIndices == 0);
        fp = size(labels_nt(strcmp(labels_nt,"target")),1);
        tn = size(labels_nt(strcmp(labels_nt,"~target")),1);
        conf_mat = [tp, fp, fn, tn];
        conf_matrices(i,:) = conf_mat;

        %calc balanced accuracy
        accuracy = (tp/(tp+fp) + tn/(tn+fn))/2;         
        %if there are no predicted targets this will give nan, so:
        if isnan(accuracy)
            accuracy = (0 + tn/(tn+fn))/2;       
        end    
        accuracies_new(i,1) = accuracy;
    end


    %sort accuracies
    [accuracies_ranked, indexOrder] = sort(accuracies_new, 'descend');

    %all genes ranked
    genenames_ranked = strings(samples,1);
    for i = 1:size(indexOrder,1)
        geneIndex = indexOrder(i);
        genesStr = string(geneNames(geneIndex));
        genenames_ranked(i) = genesStr;
    end


    %plot the accuracies overall
    figure();
    histogram(accuracies_ranked);
    title(sprintf('Accuracy for %g genes in %s', numgenes, area))
    xlabel('accuracy');
    ylabel('counts');


    %%***#1 best gene ***
    %%doesn't really make sense for cases when many genes have same top accuracy (usually
    %%numgenes > 1)
    %bestGeneIndex = indexOrder(1);
    %bestGene = string(geneNames(bestGeneIndex))
end

