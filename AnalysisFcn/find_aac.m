function rho = find_aac(Sfi,Sfj)

% Calculate amplitude amplitude correlation based on extracted spectral
% power of specific frequencies
% Sfi: freq*time
% Sfj: freq*time

N = size(Sfi,2); % number of time points
rho = zeros(size(Sfi,1),size(Sfj,1));

for i = 1:size(Sfi,1)
    for j = 1:size(Sfj,1)
        fi = Sfi(i,:);
        fj = Sfj(j,:);
        rho(i,j) = sum((fi-mean(fi)).*(fj-mean(fj)))./((N-1)*std(fi)*std(fj));
    end
end













end