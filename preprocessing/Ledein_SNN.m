function [group,coph] = Ledein_SNN(Data,system_used)
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end

% use graph based leiden algorithm
% Replace the following line by the appropriate path for Rscript
Rscript =  'D:\R\R-4.3.1\bin\Rscript.exe';
   
filefolder = 'D:\PbImpute\preprocessing';
if ~isfolder(filefolder)
    mkdir intermediateFiles
end
writetable(Data,fullfile(filefolder,'data_temporal.txt'),'Delimiter','\t','WriteRowNames',1);
% Calling R
RscriptFileName = 'D:\PbImpute\preprocessing\identify_clusters_fast.R';
%eval([' system([', '''', Rscript, RscriptFileName, '''', ' filefolder]);']);

cmd = [Rscript, ' ', RscriptFileName, ' ', filefolder];
system(cmd);

%calling R
%disp("insertint-o.R")
%RscriptFileName1 = ' ./insertinto.R ';
%eval([' system([', '''', Rscript, RscriptFileName1, '''', ' filefolder]);']);



res = 0.05:0.05:0.25; n = length(res);%3333333333333
identity_res_record = cell(1,n);
for j = 1:n
    identity = readtable(fullfile(filefolder,['identity_clustering_res',num2str(res(j)),'.txt']),'Delimiter','\t','ReadRowNames',1);
    identity_res_record{1,j} = identity.Var1;
end
[group,coph] = build_consensus(identity_res_record);

