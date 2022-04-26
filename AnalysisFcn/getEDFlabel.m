hdr = edfread('140617.edf');

newLabel = {};
for ich = 1:numel(hdr.label)
    nameInd = find(hdr.label{ich}=='-');
    newLabel{ich} = hdr.label{ich}(1:nameInd-1);
end

uniLabel = unique(newLabel,'stable')';
tchan = tabulate(newLabel);