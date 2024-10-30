function [group,coph] = partitioning(Data,M,K,numCores,system_used,accelerate,label)
if ~exist('system_used','var') || isempty(system_used)
    system_used = 'Mac';
end

if ~exist('accelerate','var') || isempty(accelerate)
    accelerate = 0;
end

if ~exist('label','var') || isempty(label)
    label = 0;
end

if label == 0
    if accelerate == 0
        poolobj = gcp('nocreate');
        if exist('numCores','var') && isempty(poolobj)
            parpool('local',numCores);
        end
        t=3; n = 20; Ds = cell(1,t);
   
        %disp(size(M'))100*264
        Ds{1,1} = pdist2(M',M','correlation'); 
        Ds{1,2} = pdist2(M',M','Spearman'); 
        Ds{1,3} = pdist2(M',M','cosin');
        %disp(size(Ds{1,1}))
 
       %222222222222222222Calculate the similarity affinity matrix.
        As = cell(1,t);
        for i = 1:t
            D = Ds{1,i}; 
            As{1,i} = exp(-D./max(D(:)));
        end
        if isempty(K)
            K  = 4;
        end
        Result = cell(1,t+1);
        for i = 1:t
            A = As{1,i};
            Re = cell(1,n);
            parfor j = 1:n
                H = symnmf_anls(A, K);
               
                Re{1,j}.H = H;
            end
            Result{1,i} = Re;
        end
        A = zeros(size(M));  A(M~=0) = 1;
        Re = cell(1,n);
        tol = 10^(-6); maxiter = 200;
        parfor i = 1:n
            disp(['INMF:',num2str(i)])
            [~,H] = INMF(M,A,K,tol,maxiter);
            Re{1,i}.H = H';
        end
        Result{1,t+1} = Re;
        %disp(Result) has four columns: columns 1 to 3 represent the affinity matrix decomposed using symNMF, and column 4 represents the INMF decomposition of the old matrix.
       % 计算
        N = size(M,2);
        Cs = cell(1,t+2);
        for flag = 1:t+1
            Re = Result{1,flag};
            label_record = cell(1,n);


         % Compute clustering labels.
            for i = 1:n
                V = Re{1,i}.H';
                label = zeros(1,N);
                for j = 1:N
                    ca = find(V(:,j) == max(V(:,j)));%The index of the maximum value in column j, which corresponds to the predicted clustering label.
                    %For example, if the maximum value in the first column is in the third row, then the cell belongs to cluster 3.
                    label(j) = ca(1);
                end
                label_record{1,i} = label;
            end
      

            % compute consensus matrix
            consensus = zeros(N);
            for i = 1:n
                label = label_record{1,i}; C = zeros(N);
                for j = 1:N
                    for k = j:N
                        if label(j) == label(k)
                            C(j,k) = C(j,k)+1;
                            C(k,j) = C(j,k);
                           % Connectivity matrix C
                        end
                    end
                end
                consensus = consensus + C;
            end
            Cs{1,flag} = consensus/n;
        end
   % The last one is the sum of the first three.
        C = zeros(N);
        for i = 1:t
            C = C + Cs{1,i};
        end
        Cs{1,t+2} = C/t;
  




        % compute the coph of consensus matrix()
        clusts = cell(1,t+2); cophs = zeros(1,t+2);
        for i = 1:t+2
            [~,clusts{1,i},~,cophs(i)] = nmforderconsensus0(Cs{1,i},K);
            %
        end
        cophid = cophs(t+1:t+2);
        if abs(max(cophid)-min(cophid)) > 0.1
            id = find(cophid == max(cophid));
            group = clusts{1,t+id}; coph = cophs(t+id);
        else
            C = (Cs{1,t+1}+Cs{1,t+2})/2;
            [~,group,~,coph] = nmforderconsensus0(C,K);
        end
        % delete(gcp('nocreate'))
        disp("Clusters of the clustering")
        disp(group)
        disp("Clustering performance")
        disp("coph:"+coph)
        % testtsne();
    else
        [group,coph] = Ledein_SNN(Data,system_used);
        disp("Clusters of the clustering")
         
fileID = fopen('group.txt', 'w');
fprintf(fileID, '%d\n', group);
fclose(fileID);

        %testtsne();
    end
else
    group = label; coph = 1;
end

