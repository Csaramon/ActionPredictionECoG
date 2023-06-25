function hh = goodbar(barH,barE,pMatrix,barX)


if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

if nargin < 2
    barE = [];
    pMatrix = [];
    barX = 1:numel(barH);
end

if nargin < 3
    pMatrix = [];
    barX = 1:numel(barH);
end
if nargin < 4
    barX = 1:numel(barH);
end

pMatrix = round(pMatrix*1000)/1000;
pMatrix(pMatrix==0) = 0.001;

if size(barH,1) > size(barH,1)
    barH = barH';
end
% parameter for error bar tick
tickLength = 0.1*min(pdist(barX'));
tickColor = [0,0,0]; % black
tickWidth = 1.5;

hold(gca,'on')
for ib = 1:numel(barH)
    hh(ib) = bar(barX(ib),barH(ib),'EdgeAlpha',1,'LineWidth',tickWidth,'FaceColor',[1 1 1]);
    if ~isempty(barE)
        if barH(ib) >= 0
        line([barX(ib),barX(ib)],[barH(ib),barH(ib)+barE(ib)],'color',tickColor,'linewidth',tickWidth)
         line([barX(ib)-tickLength/2,barX(ib)+tickLength/2],[barH(ib)+barE(ib),barH(ib)+barE(ib)],'color',tickColor,'linewidth',tickWidth)
        elseif barH(ib) < 0
            line([barX(ib),barX(ib)],[barH(ib),barH(ib)-barE(ib)],'color',tickColor,'linewidth',tickWidth)
             line([barX(ib)-tickLength/2,barX(ib)+tickLength/2],[barH(ib)-barE(ib),barH(ib)-barE(ib)],'color',tickColor,'linewidth',tickWidth)
        end
       
    end
    
end

addH = zeros(size(barH));

if ~isempty(pMatrix)
    if all(size(pMatrix) == [numel(barH),numel(barH)])
        for i = 1:size(pMatrix,1)
            for j = 1:size(pMatrix,2)
                if ~isnan(pMatrix(i,j))
                    if i==j
                         text(barX(i)-2*tickLength,1.05*barH(i)+max(barE),['p = ' num2str(pMatrix(i,j))])
                    else
%                         line([barX(i),barX(j)],[1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j)), ...
%                             1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j))],'color',tickColor,'linewidth',tickWidth)
%                         line([barX(i),barX(i)],[1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j))-0.015, ...
%                             1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j))],'color',tickColor,'linewidth',tickWidth)
%                         line([barX(j),barX(j)],[1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j))-0.015, ...
%                             1.04*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j))],'color',tickColor,'linewidth',tickWidth)
                        text(mean([barX(i),barX(j)]-3*tickLength),1.075*max([barH(i),barH(j)])+max(barE)+max(addH(i),addH(j)),['p = ' num2str(pMatrix(i,j))])
                        addH([i,j]) = addH([i,j])+0.04*max(barH);
                    end
                end
            end      
        end
        
        
        
    else
        error(message('Incorrect size of p value matrix.'));
    end
end

% ylim([0.85*min(barH),1.15*max(barH)])

end

