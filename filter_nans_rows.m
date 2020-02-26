clear all;
close all;
load("C:\Users\elobe\Dropbox\1_HONOURS\code\data\AllenGeneDataset_19419_Ben.mat", '-mat')

genes = GeneExpData.combZ.energy;

[rows, cols] = size(genes);

%nans in each row
nans_in_rows = zeros(1,rows);
for i = 1:rows
row = genes(i,:);
[x, nans_row] = size(row(isnan(row)));
nans_in_rows(i) = nans_row;
end

%now see how many we would be deleting in the cases where we
%delete anything with more than 1%, 5%, 10%, 15% nans

perc_vals = [1,5,10,15,20,25,30,35,40,50];
rows_deleted = [];
perc_rows_deleted = [];
for i = 1:length(perc_vals)
    [x,r_d] = size(nans_in_rows(nans_in_rows > perc_vals(i)));
    p_r_d = r_d/rows*100;
    rows_deleted = [rows_deleted r_d];
    perc_rows_deleted = [perc_rows_deleted p_r_d];
end

%view rows_deleted and perc_rows_deleted
