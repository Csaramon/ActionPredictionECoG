allColors = distinguishable_colors(numel(Labels));
allColors = round(allColors.*255);

    fid=fopen(['AnatomyMacroLUT.txt'],'w+','native','UTF-8');
    for i = 1:numel(Labels)
        fprintf(fid,'%-5d %-30s %-5d %-5d %-5d %-5d\n',...
            i,Labels{i},allColors(i,1),allColors(i,2),allColors(i,3),0);
    end
    fclose(fid);
