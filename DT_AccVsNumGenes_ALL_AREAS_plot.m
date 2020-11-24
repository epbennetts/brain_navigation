close all
figure();
% %set(0,'DefaultAxesColorOrder',hsv(numAreas));
% %which areas to plot: all -> 1:numAreas, particular area --> index in areaNames
% areas_bigger = [1,2,3,5,7,8,9,10,11]; %those with >10 data points
% areas_plotting = areas_bigger;
% numAreas_plotting = size(areas_plotting);
% 
% %make area labels match the areas being plotted
% areaNames_sel = areaNames(areas_plotting);


%for a = 1:numAreas_plotting
allIndices = 1:numAreas;
areasBiggest_indices = [1,7,9,11];
Isocortex_ind = [1];
areasIncrease_indices = [7,9,11]; %areas that go to more than 1 gene

%set which areas to display
areaIndices = Isocortex_ind;
areaNames_disp = areaNames(areaIndices);

%which stopping criterion? 
stoppingCritSetting = 'no';

%show legend?
showLegend = 'no';

for a = areaIndices
    area = areaNames(a);
    bestAccs_temp = bestAccs_ALL(:,a);
    bestAccs_stdevs_temp = bestAccs_stdevs_ALL(:,a);
    
%     %SIMPLEST STOPPING CRITERION
%     for i = 2:size(bestAccs_temp)
%        %for noise = 0 case: if saturates,  make the rest of array NaN 
%        if (bestAccs_temp(i) <= bestAccs_temp(i-1))
%            bestAccs_temp(i:size(bestAccs_temp)) = NaN;
%        end
%     end
    
%     %PROVISIONAL STOPPING CRIT FOR SPECIFIC PLOT (for now, do properly later)
%     if a == 7 || a == 9 || a == 11
%         bestAccs_temp(3:size(bestAccs_temp)) = NaN;
%     else
%         bestAccs_temp(2:size(bestAccs_temp)) = NaN;
%     end
    
%PROPER STOPPING CRIT
if (strcmp(stoppingCritSetting,'yes'))
    for i = 2:size(bestAccs_temp)
        %for noise = 0 case: if saturates,  make the rest of array NaN
        if (bestAccs_temp(i)-bestAccs_stdevs_temp(i)  <= bestAccs_temp(i-1)+bestAccs_stdevs_temp(i-1))
            bestAccs_temp(i:size(bestAccs_temp)) = NaN;
        end
    end
end
    
    
    
    %Accuracies vs numgenes    
    %Plot accuracy increase:
    numgenes_array = 1:maxNumGenesInDT;
    p1 = errorbar(numgenes_array, bestAccs_temp.*100, bestAccs_stdevs_temp.*100,'.-');
    %p1 = plot(numgenes_array, bestAccs_temp, '.-', 'LineWidth',1);
    p1.Color(4) = 0.4;
    %legend(areaNames, 'Location', 'bestOutside')
    if strcmp(showLegend, 'yes')
        legend(areaNames_disp, 'Location', 'bestOutside')
    end
    %title(sprintf('Balanced Accuracy vs num genes for %d noise (%d samples, %d genes, %d iters)', noiseStDev, sizeSampleSubset, maxNumGenesInDT, numNoiseIterations));
    %s (samples: %d, iters: %d, folds: %d, numGenes: %d)', area, sizeSampleSubset, numNoiseIterations, numFolds, numGenes_temp))
    xlabel('Sensors')
    ylabel('Balanced Accuracy (%)')
    set(gca,'xtick', 0:maxNumGenesInDT)
    xlim([0.5 10.5])
    ylim([70 101])
    grid on;
    hold on;
end

%legend('areas')