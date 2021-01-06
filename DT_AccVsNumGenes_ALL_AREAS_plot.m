%close all
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
areasIncrease_indices = [11,7,9]; %areas that go to more than 1 gene

%set which areas to display
areaIndices = areasIncrease_indices;
%areaIndices = Isocortex_ind;
areaNames_disp = areaNames(areaIndices);
numAreas_disp = size(areaIndices,2);

%empty arrays
bestAccs_stopping = NaN(numAreas,1); %<------------- make array of accuracies at each stopping crit
errors_stopping = NaN(numAreas,1);

% %colours (Isocortex)
% colours = {'009FAC'};
%colours (for the 3 areas)
colours = {'FF909F', 'FFA6FF', 'FFB3D9'};

%which stopping criterion?
stoppingCritSetting = 'no';

%show legend?
showLegend = 'yes';

%plot bar plot as well?
showBarPlot = 'no';

%plot acc vs numgenes for multiple areas
counter = 1;
for a = areaIndices
    area = areaNames(a);
    bestAccs_temp = bestAccs_ALL(:,a);
    bestAccs_stdevs_temp = bestAccs_stdevs_ALL(:,a);
    colour = colours{1,counter};
    counter = counter+1;
    
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
            if (bestAccs_temp(i)  <= bestAccs_temp(i-1)+bestAccs_stdevs_temp(i-1))
                %save accuracies for each area at stopping crit
                bestAccs_stopping(a) = bestAccs_temp(i-1);
                errors_stopping(a) = bestAccs_stdevs_temp(i-1);
                
                %don't plot following stopping crit
                %bestAccs_temp(i:size(bestAccs_temp)) = NaN;
            end
        end
    end
    
    if isnan(bestAccs_stopping(a))
        disp('using acc from gene number 10')
        bestAccs_stopping(a) = bestAccs_temp(10);
        errors_stopping(a) = bestAccs_stdevs_temp(10);
    end
    
    
    
    %Accuracies vs numgenes
    %Plot accuracy increase:
    numgenes_array = 1:maxNumGenesInDT;
    p1 = errorbar(numgenes_array, bestAccs_temp.*100, bestAccs_stdevs_temp.*100,'.-');
    %p1 = plot(numgenes_array, bestAccs_temp, '.-', 'LineWidth',1);
    %p1.Color(4) = 0.4;
    %p1.Color = rgbconv(colour);
    %p1.Color = rgbconv('1F9D5A');%'009FAC');%'248A5E'); %Isocortex
    p1.LineWidth = 1;
    p1.MarkerSize = 12;
    set(gca, 'FontSize', 13);
    %legend(areaNames, 'Location', 'bestOutside')
    if strcmp(showLegend, 'yes')
        legend(areaNames_disp, 'Location', [0.65, 0.16, 0.2, 0.2]) %all
        %legend(areaNames_disp, 'Location', [0.65, 0.25, 0.2, 0.1]) % Isocortex
    end
    
    %title(sprintf('Balanced Accuracy vs num genes for %d noise (%d samples, %d genes, %d iters)', noiseStDev, sizeSampleSubset, maxNumGenesInDT, numNoiseIterations));
    %s (samples: %d, iters: %d, folds: %d, numGenes: %d)', area, sizeSampleSubset, numNoiseIterations, numFolds, numGenes_temp))
    xlabel('N: Num genes')
    ylabel('BA: Balanced Accuracy [%]')
    set(gca,'xtick', 0:maxNumGenesInDT)
    xlim([0.5 10.5])
    ylim([65 100])
    grid on;
    hold on;
end
hold off;

%convert to sight subset
bestAccs_stopping = bestAccs_stopping(areaIndices);
errors_stopping = errors_stopping(areaIndices);


%legend('areas')

%plot bar plot of accuracies in different areas
if strcmp(showBarPlot, 'yes')
    figure();
    bestAccs_temp2 = bestAccs_ALL(10,:).*100;
    bestAccs_stdevs_temp2 = bestAccs_stdevs_ALL(10,:).*100;
    
%     X = categorical(areaNames_disp);
%     Y = bestAccs_temp2(1,areaIndices)';
%     errors = bestAccs_stdevs_temp2(1,areaIndices);
     X = categorical(areaNames_disp);
     Y = bestAccs_stopping.*100; 
     errors = errors_stopping.*100;
    
    
    %bar plot
    bar(X,Y);
    %title(sprintf('Accuracy vs area'));
    %xlabel('Brain Area')
    ylabel('Balanced Accuracy (%)')
    hold on;
    
    %error bars
    %errorbar(Y,errors);    %<-----------------------------
    %     er.Linestyle = 'none';
    
    
    %set(gca,'xtick', 0:numgenes)
    grid on;
    hold on;
    
end

