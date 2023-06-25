% cluster based permutation test subfcn
    function [cStat,cStatAll,L,cSTATS] = calculateClusterStat(dataA,dataB,pvalue,metric)
        
        
        [Hi,~,~,cSTATS] = ttest(dataA,dataB,'alpha',pvalue,'tail','both');%,'tail','right'
    
        tpic = Hi.*cSTATS.tstat;
        bipic = (tpic~=0);       %bianarization
        L = bwconncomp(bipic,8);
        switch metric
            case 'maxsum'
                statistic_tMap= regionprops(L,tpic,'PixelValues');
                cStatAll = zeros(size(statistic_tMap));
                
                for is = 1:size(statistic_tMap,1)
                    
                    cStatAll(is) = sum(statistic_tMap(is).PixelValues);
                    
                    
                end
                cStat = max(abs(cStatAll));
                
                
            case 'maxsize'
                cluster_area = regionprops(L,'Area');
                cStat = max(cat(1,cluster_area.Area));
                
                
        end
        
        % no significant cluster
        if isempty(cStatAll)
            cStat = 0;
            cStatAll = 0;
        end
        
    end