function id = gene_selection(M0,iniData,flag,system_used,low_mu,high_mu,low_F)
if ~exist('flag','var') || isempty(flag)
    flag = 0; % select gene based on Fano factor and the mean expression
end
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac'; % select gene based on Fano factor and the mean expression
end
    Rscript = 'D:\R\R-4.3.1\bin\Rscript.exe'; % for 64-bit windows
if ~exist('low_mu','var') || isempty(low_mu)
  
    low_mu = 0.01;
end
if ~exist('high_mu','var') || isempty(high_mu)
     %high_mu = 3.5;
     high_mu = 3.5;
end
if ~exist('low_F','var') || isempty(low_F)
     %low_F = 0.01;
     low_F = 0.01;
end

% disp('selecting genes:');
if flag == 0
    disp('fano')
    id = HVGs(M0,low_mu,high_mu,low_F,iniData);%fanolow_mu = 0.01;high_mu = 3.5;  low_F = 0.5;
    selected_genes_matrix = iniData(id, :);
    disp(size(selected_genes_matrix));
    
    

% % Specify the file path and filename to save
selected_genes_file_path = 'D:\PbImpute\preprocessing\selected_genes_expression_matrix.csv';
% 
% % Use the writetable function to save the selected gene expression matrix as a .txt file
writetable(selected_genes_matrix, selected_genes_file_path);
% 
% 



elseif flag == 1
    filefolder = pwd;%Get the current directory
    T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
    writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
    % Calling R's GiniIndex
    eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
    id = importdata('Gini_ID.txt');

else
    id1 = HVGs(M0,low_mu,high_mu,low_F);
    T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
    writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
    filefolder = pwd;
    eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
    id2 = importdata('Gini_ID.txt');
    id = union(id1,id2);
end


