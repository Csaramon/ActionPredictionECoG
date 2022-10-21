function BF = BCT(Dat,Hypo)

% function input:
% Dat: m*n matrix containing m observations and n variables
% Hypo: Hypothesis to test in the format of 'X1_with_X2 > X1_with_X3 > 0'
%                in which X1...Xn represent the Nth variable, use semicolon
%                to separate multiple hypotheses
%                
if nargin < 2
    Hypo = [];
end
if nargin < 1
    warning('No input data! A result based on random matrix will be generated.')
    Dat = rand(100,5);
end

P = mfilename('fullpath');
mfilepath = strrep(P,mfilename,'');

% Create a temporary matrix for R script
save([mfilepath 'tmp.mat'],'Dat','-v6')

if strcmpi(computer,'PCWIN64')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox')
    addpath('C:\Users\qin2\Documents\MATLAB\toolbox\fieldtrip-20210418')
    basePath = 'C:\Users\qin2\Documents\ActionPredictionECoG\';
elseif strcmpi(computer,'MACI64')
[~,cmdResult] = unix(['. ~/.bashrc;. ~/.zshrc;' ...
    'export Hypo="' Hypo '";'...
    'Rscript test.r ' mfilepath 'tmp.mat $Hypo']);
elseif strcmpi(computer,'GLNXA64')
    addpath('/data00/Chaoyi/ActionPredictionECoG/')
    addpath('/data00/Chaoyi/toolbox/fieldtrip-20210418/')
    basePath = '/data00/Chaoyi/ActionPredictionECoG/';
end
delete([mfilepath 'tmp.mat'])

TestInd = strfind(cmdResult,'Bayesian hypothesis test');
BF.Output = cmdResult(TestInd(1):end);

% Format a bit of the output for future use
aa = strsplit(BF.Output,'\n')';
PrStartInd = find(strcmp('Posterior probabilities:',aa),1);
PrEndInd = find(strcmp('Bayesian hypothesis test',aa),1,'last');
if PrEndInd > PrStartInd
    tmpCell = cellfun(@(x) strsplit(x,' '),aa(PrStartInd+1:PrEndInd-1),'UniformOutput',false);
    BF.Pr = vertcat(tmpCell{:});
else
    tmpCell = cellfun(@(x) strsplit(x,' '),aa(PrStartInd+1:end-1),'UniformOutput',false);
    BF.Pr = vertcat(tmpCell{:});
    return
end
PrhStartInd = find(strcmp('Posterior probabilities:',aa),1,'last');
EvidenceInd = find(strcmp('Evidence matrix (Bayes factors):',aa),1);
tmpCell = cellfun(@(x) strsplit(x,' '),aa(PrhStartInd+1:EvidenceInd-1),'UniformOutput',false);
BF.Prh = vertcat(tmpCell{:});

SpecificationInd = find(strcmp('Specification table:',aa),1);
tmpCell = cellfun(@(x) strsplit(x,' '),aa(EvidenceInd+1:SpecificationInd-1),'UniformOutput',false);
BF.Evidence = vertcat(tmpCell{:});

HypothesesInd = find(strcmp('Hypotheses:',aa),1);
tmpCell = cellfun(@(x) strsplit(x,' '),aa(SpecificationInd+1:HypothesesInd-1),'UniformOutput',false);
BF.Specification = vertcat(tmpCell{:});
BF.Hypotheses = aa(HypothesesInd:end-1);

end

