function [hf,ha,h1,h2,rho,pvalue] = plotCorrelation(data1,data2)



[rho,pvalue] = corr(data1,data2);   

p = polyfit(data1,data2,1);
ydata = polyval(p,data1);
hf = figure('Name','Correlation');
ha = axes('Parent',hf);


h1 = plot(ha,data1,data2,'.','markersize',10);
hold on
h2 = plot(ha,data1,ydata);


end


