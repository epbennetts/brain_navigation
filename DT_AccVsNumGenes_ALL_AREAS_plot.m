

figure();
%set(0,'DefaultAxesColorOrder',hsv(numAreas));
for a = 1:numAreas
    area = areaNames(a);
    bestAccs_temp = bestAccs_ALL(:,a);
    bestAccs_stdevs_temp = bestAccs_stdevs_ALL(:,a);
    
    
    for i = 2:size(bestAccs_temp)
       %for noise = 0 case: if saturates,  make the rest of array NaN 
       if noiseStDev==0 && (bestAccs_temp(i) == bestAccs_temp(i-1))
           bestAccs_temp(i:size(bestAccs_temp)) = NaN;
       end
    end
    
    %Accuracies vs numgenes    
    %Plot accuracy increase:
    numgenes_array = 1:maxNumGenesInDT;
    p1 = errorbar(numgenes_array, bestAccs_temp, bestAccs_stdevs_temp,'.-');
    %p1 = plot(numgenes_array, bestAccs_temp, '.-', 'LineWidth',1);
    p1.Color(4) = 0.4;
    legend(areaNames, 'Location', 'bestOutside')
    title(sprintf('Balanced Accuracy vs num genes in all areas (noise: %d)', noiseStDev));
    xlabel('Num of genes used in Decision Tree')
    ylabel('Balanced Accuracy (%)')
    set(gca,'xtick', 0:maxNumGenesInDT)
    xlim([0.5 10.5])
    ylim([0.7 1.01])
    grid on;
    hold on;
end

%legend('areas')