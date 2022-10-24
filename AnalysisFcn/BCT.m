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

% set up and run R in terminal
if strcmpi(computer,'PCWIN64')
    if exist('C:\Program Files\R\R-4.2.1\bin\Rscript.exe','file')
        setenv('PATH', 'C:\Program Files\R\R-4.2.1\bin\') % default path for R, change it for your own path
    else
        error('Unable to find R in default location, please specify your own path to R')
    end
    [~,cmdResult] = system(['Rscript ' mfilepath 'BCT.R ' mfilepath 'tmp.mat "' Hypo '"']);
elseif strcmpi(computer,'MACI64')
    [~,cmdResult] = unix(['. ~/.bashrc;. ~/.profile;. ~/.zshrc;' ...
        'export Hypo="' Hypo '";'...
        'Rscript ' mfilepath 'BCT.R ' mfilepath 'tmp.mat $Hypo']);
elseif strcmpi(computer,'GLNXA64')
    if isempty(Hypo)
        [~,cmdResult] = system(['Rscript BCT.R ' mfilepath 'tmp.mat']);
    else
        [~,cmdResult] = system(['export Hypo="' Hypo '";'...
            'Rscript ' mfilepath 'BCT.R ' mfilepath 'tmp.mat "' Hypo '"']);
    end
end
% delete tmp data file
delete([mfilepath 'tmp.mat'])

% get the original output from R
TestInd = strfind(cmdResult,'Bayesian hypothesis test');
BF.Output = cmdResult(TestInd(1):end);

% Format a bit of the output for other purpose
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

