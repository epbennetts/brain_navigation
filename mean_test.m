% % % % % 
% % % % % x = [5 2 7 5 8];
% % % % % X = sort(x, 'descend');
% % % % % [Y, Z] = sort(x, 'descend');
% % % % % 
% % % % % 
% % % % 
% % % % [genes, targetIndices, classes, geneNames] = filter_nans;
% % % % target = genes()
% % % 
% % % 
% % % county = 0;
% % % for i = [1,3,4,6,8]
% % %     county = county + 1;
% % %    figure(i); 
% % % end
% % % 
% % % % close all;
% % % % a = [1,3,4,6,8];
% % % % county = 0;
% % % % for i = a
% % % %    county = county + 1;
% % % %    figure(i);  
% % % % end
% % 
% % clear vars;
% % a = "blah";
% % b = "blah";
% % strcmp(a,b) 
% % 
% % if (strcmp(a,b))
% %     disp("woop woop");
% % end
% % 
% 
% 
% X = [0,0,0;
%     1,1,1;
%     2,2,2];
% a = mean(X,1)
% b = mean(X,2)



numgenes_array = 1:numgenes;
figure();
plot(numgenes_array, top_accuracy,'.-b');
title(sprintf('Accuracy vs # genes in %s', area));
xlabel('Num genes used in DTs')
ylabel('Accuracy')
set(gca,'xtick',0:numgenes)
grid on;
