function [cluster,coph] = build_consensus(identity_res_record)
n = length(identity_res_record); Cs = cell(1,n);
for i = 1:n
    identity = identity_res_record{1,i};%
    N = length(identity);
    % compute consensus matrix
    C = zeros(N);
    for j = 1:N
        for k = j:N
            if isequal(identity(j),identity(k))
                C(j,k) = C(j,k)+1;
                C(k,j) = C(j,k);
            end
        end
    end
    Cs{1,i} = C;
end
C = zeros(N);
for j = 1:n
    C = C+Cs{1,j};
end
C = C/n;
% determine k
nPC = 30;
[~,S,~] = rsvd(C,nPC); s = diag(S); s = [s;zeros(nPC-length(s),1)];
% compute ratio
% cs = cumsum(s)/sum(s);
tol=0.15;
disp("666666666666")
disp((s/sum(s)))
K=3%test K=3
%K = find((s/sum(s) < tol) == 1, 1 )-1;

% figure; plot(1:length(cs),cs,'LineWidth',1.5)
% ax=gca; ax.LineWidth=1;
% set(gca,'Xtick',1:length(cs))
% label = cell(1,length(cs));
% for i = 1:length(cs)
%     label{i} = num2str(i);
% end
% set(gca,'xticklabel',label);
% hold on;
% th = 0:pi/50:2*pi;
% xunit = 0.1 * cos(th) + K;
% yunit = 0.1 * sin(th) + cs(K);
% plot(xunit, yunit);
% hold off
[~,cluster,~,coph] = nmforderconsensus0(C,K);
disp(coph)