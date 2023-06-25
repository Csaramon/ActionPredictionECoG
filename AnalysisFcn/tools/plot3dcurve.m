function hl = plot3dcurve(p0,p1,colorvec,linewidth)
if nargin < 4
    linewidth = 2;
end
if nargin < 3
    colorvec = [1 0 0];
end
t = 0:0.01:1;
n = 1;
H=3;
for j = t
    a(n)=p0(1)+(p1(1)-p0(1))*j;
    b(n) = p0(2)+(p1(2)-p0(2))*j;
    h(n) = p0(3)+(p1(3)-p0(3))*j+H*sin(pi*j);
    n= n + 1;
end
hl = line(a,b,h,'Color',colorvec,'LineWidth',linewidth);