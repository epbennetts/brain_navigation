    %TO DO
    %later: have clean genes matrix saved somewhere (instead of processing every time)
    %generalise to a function that uses n num of genes and just goes through
    %everything without plotting
    %generalise further to any area

%     numgenes = 1;
%     samples = 19114;
%     area = "Isocortex";
%     prev_best_genes = [1237, 7, 48];


function [indexOrder, accuracies_ranked, genenames_ranked, thresholds_all] = DT_classification_multiple(prev_best_genes, samples, area)
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


    numgenes = length(prev_best_genes);

    % COULD BE EASIER?:
    % isTarget
    % targetIndices = find(isTarget);
    % notTargetIndices = find(~isTarget);

    [genes, targetIndices, target, nonTarget, classes, geneNames] = filter_nans(area);
    [rows, cols] = size(genes);

    %num of trees we want to plot (not in order)
    %plotting = 1:10;
    %plotting = [];

    %set vars and empty arrays
    accuracies_new = zeros(samples,1);
    f1_score = zeros(samples,1);
    conf_matrices = zeros(samples,4);
    thresholds_all = NaN(samples, numgenes);
    %genes_used = cell(samples, 2);

    % Previously selected genes as baseline:
    prevGeneData = genes(:,prev_best_genes);

    %loop, make trees, classify and evaluate
    for i = 1:samples

        % Set gene data for this iteration:
        gene_combo = [prevGeneData, genes(:,i)];

        %train
        tree = fitctree(gene_combo, classes, 'MaxNumSplits', numgenes);

        %check first predictor (should be best gene: x1)
        %doesn't exactly work as doesn't always show second split?
        %gene_used = tree.CutPredictor;
        %genes_used{i,1} = gene_used{1};
        %genes_used{i,2} = gene_used{2};

        % Store all thresholds for later:
        thresholds_raw = tree.CutPoint; %threshold with NaNs
        thresholds = thresholds_raw(~isnan(thresholds_raw));
        threshold_pan = nan(1,numgenes);
        threshold_pan(1:length(thresholds)) = thresholds;
        thresholds_all(i,:) = threshold_pan;
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

        % MAYBE CLEARER WITH BINARY FORMULATION:
        % tp = sum(labels==1 & isTarget==1);
        % fn = sum(labels==1 & isTarget==0);

        % labels_t is the set of predictions for areas that are targets
        labels_t = labels(targetIndices == 1);
        tp = sum(strcmp(labels_t,'target')); % much faster if binary labels
        fn = sum(strcmp(labels_t,'~target');
        %non-targets
        labels_nt = labels(targetIndices == 0);
        fp = sum(strcmp(labels_nt,'target'));
        tn = sum(strcmp(labels_nt,'~target'));
        conf_matrices(i,:) = [tp, fp, fn, tn];

        %calc balanced accuracy
        accuracies_new(i) = (tp/(tp+fp) + tn/(tn+fn))/2;
        prec = tp/(tp+fp);
        recall = tp/(tp+fn);
        f1_score(i) = 2*prec*recall/(prec+recall);
        %if there are no predicted targets this will give nan, so:
        % if isnan(accuracy)
        %     accuracy = (0 + tn/(tn+fn))/2;
        % end
    end


    %sort accuracies
    [accuracies_ranked, indexOrder] = sort(accuracies_new, 'descend');

    % indexOrder ~ [4,2,3,1]
    % X = [5,6,2,1]; X(indexOrder) => [1,6,2,5]

    % geneNamesRanked = genenames{indexOrder};

    % %all genes ranked
    % genenames_ranked = strings(samples,1);
    % for i = 1:size(indexOrder,1)
    %     geneIndex = indexOrder(i);
    %     genesStr = string(geneNames(geneIndex));
    %     genenames_ranked(i) = genesStr;
    % end


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
