%TO DO 
%put all cols into DT
%later: have clean genes matrix saved somewhere (instead of processing every time)

[genes, targetIndices] = filter_nans;
[rows, cols] = size(genes);


tree = fitctree(genes(:,1),targetIndices, 'Crossval', 'on');
view(tree.Trained{1},'Mode','graph')

    
  