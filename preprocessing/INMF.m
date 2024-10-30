function [U,V] = INMF(X,A,k,tol,maxiter,varargin)
% used for the following model:
% min (||A.*(X-UV)||_F)^2
% <Inputs>:
% X: the real full data matrix
% A: the observed projected data matix, with elsements equals 0 represents
% missing data index, 1 represents observed data
% k: the low rank
% tol: the stop tolerennce
% maxiter: the maximum iterations
% <Outputs>:
% U: the loading matrix
% V: the coefficent matrix
%% main function
[m,n] = size(X);
% initialization
U = rand(m,k); V = rand(k,n);
% Read optional parameters
if (rem(length(varargin),2) == 1)
    error('Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'U_INIT',      U = varargin{i+1};
            case 'V_INIT',      V = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end

for i = 1:maxiter
    U = (U.*((A.*X)*V'))./((A.*(U*V)*V')+eps); % update U
    V = V.*(U'*(A.*X))./(U'*(A.*(U*V))+eps);
    stop_value = norm(A.*(X-U*V),'fro')^2; 
    if stop_value < tol
        break;
    end
end

    
    
