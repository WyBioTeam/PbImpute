%% %% This code defines a function named HVGs, which selects high variable genes (HVGs) based on the mean gene expression and Fano factor.
function id = HVGs(M,low_mu,high_mu,low_F,iniData)
if ~exist('low_mu','var') || isempty(low_mu)
    low_mu = 0; % select gene above this average expression level
end
if ~exist('high_mu','var') || isempty(high_mu)
    high_mu = 10; % select gene below this average expression level
end
if ~exist('low_F','var') || isempty(low_F)
    low_F = 0; % select gene above this Fano factor
end
% M is the raw count data matrix
% id is the set of select gene position

%% select HVGs: calculates the average expression and Fano factor for each gene, places these genes into bins,
% and then calculates a z-score for Fano factor within each bin.
mu = log(1+mean(10.^M-1,2));%Perform averaging calculations, where the original data is logarithmic and is converted back to the original scale before averaging, and then converted back to logarithmic. Parameter 2 indicates that the operation is performed across each row (each gene class).

F = log(var(10.^M-1,0,2)./mean(10.^M-1,2));
mu(isnan(mu)) = 0;F(isnan(F)) = 0;
mu(isinf(mu)) = 0;F(isinf(F)) = 0;
%% Divide the average expression values, mu, into 20 discrete bins, and return the bin labels Y and boundaries E.
[Y,E] = discretize(mu,20);%Binning process: Y represents the bin number to which each value in the mu data belongs, and E contains 21 range boundaries
% E represents the boundaries as follows: minimum value ---> (maximum - minimum) / 20 ---> 2 * (maximum - minimum) / 20 ---> 3 * (maximum - minimum) / 20 ... ---> 19 * (maximum - minimum) / 20 ---> maximum value.
idx = setdiff(1:20,unique(Y));
if ~isempty(idx)
     %If the index does not cover all bins from 1 to 20, then re-bin the data, for example, if index 5 is missing (to prevent clustering without genes).
    E(idx+1) = [];
    Y = discretize(mu,E);
end
%% Calculate the mean and standard deviation of the Fano factor within each bin.
mean_y = grpstats(F,Y,"mean");
sd_y = grpstats(F,Y,"std");
F_scaled = (F - mean_y(Y))./sd_y(Y);% Normalize each bin.
F_scaled(isnan(F_scaled)) = 0;
id = find((mu > low_mu) & (mu < high_mu) & F_scaled > low_F);
disp("ggggggggggggggg")
disp(length(id))
%id = sortrows([id, F_scaled(id)], 2, 'descend');  % Sort genes based on Fano factor
%id = id(1:2000, 1);  % Select top 2000 genes
%id = sort(id); % 
%disp(id)
%length(id)
figure
scatter(mu,F_scaled,'k.')
xlabel('Average expression');
ylabel('Fano factor')

%bulkdata = readtable('bulkFetalBrain.csv'); 
%selected_bulk = bulkdata(id, :);

%writetable(selected_bulk, 'bulkFetalBrain_2000.csv'); 


selected_M = iniData(id, :); 


