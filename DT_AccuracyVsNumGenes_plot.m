%Uses the .mat file or workspace obtained from cluster to plot Accuracy vs numGenesC 


% ACCURACIES PLOT
DT_plot_accuracies(maxNumGenesInDT, top_accuracies, area)

%PLOT GENE EXPRESSION & THRESHOLDS? 
% later in separate function:
% %Plot thresholds
% numplots = 5;
% range = 'top';
% %work on this functionethi
% %DT_plot_multiple(genes, geneNames, numgenes, best_genes, isTarget, thresholds_all, samples, indexOrder, numplots, range)

%end of program
% WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
% Audio = audioplayer(WarnWave, 22050);
% play(Audio);