function [FDR,Q] = myfdr(p,alpha)

% calulate fdr using Benjamini and Hochberg's method

if nargin < 2
    alpha = 0.05;
end


[sortP,pIndex] = sort(p);
FDR = zeros(size(p));
Q = zeros(size(p));

for irank = 1:numel(sortP)


    FDR(pIndex(irank)) = alpha*irank/numel(sortP);
    Q(pIndex(irank)) = sortP(irank)*numel(sortP)/irank;


end