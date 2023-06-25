function dottedErrorbar(meanVal,seVal,x)

if nargin < 3
    x = 1:numel(meanVal);
end
figure()
xlim([min(x)-1,max(x)+1])
hold on
plot(x,meanVal,'ko','markersize',10,'linewidth',1.5)
plot(x,meanVal,'k','linewidth',1.5)

for ix = 1:numel(x)
plot([x(ix),x(ix)],[meanVal(ix)-seVal(ix),meanVal(ix)+seVal(ix)],'k','linewidth',1.5)

plot([x(ix)-0.1,x(ix)+0.1],[meanVal(ix)-seVal(ix),meanVal(ix)-seVal(ix)],'k','linewidth',1.5)
plot([x(ix)-0.1,x(ix)+0.1],[meanVal(ix)+seVal(ix),meanVal(ix)+seVal(ix)],'k','linewidth',1.5)


end





end