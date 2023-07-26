function data = find_doublets(data,dt,flag,N)
%
% Exclude/include doublets in a ripple raster plot
% Input:
% data = time * trials ripple raster
% dt = doublet max interval in samples
% flag = 0 - exclude doublets / 1 = include only doublets
% N = number of events (2=doublets, 3=triplets, 4=quadruplets)
% Author: Itzik Norman 2019

if ~exist('dt','var'), dt = 150; end
if ~exist('flag','var'), flag = 0; end
if ~exist('N','var'), N = 2; end

switch N
    case 2
        ind2 = find_doublet_indices(data,dt); % doublets
        ind = ind2;
    case 3
        ind2 = find_doublet_indices(data,dt);
        datatmp = zeros(size(data)); datatmp(ind2) = 1;
        ind3 = find_doublet_indices(datatmp,dt); % triplets
        ind = ind3;
    case 4
        ind2 = find_doublet_indices(data,dt);
        datatmp = zeros(size(data)); datatmp(ind2) = 1;
        ind3 = find_doublet_indices(datatmp,dt);
        datatmp = zeros(size(data)); datatmp(ind3) = 1;
        ind4 = find_doublet_indices(datatmp,dt); % quadruplets
        ind = ind4;
    otherwise
        error('Wrong number of events (check the argument ''N'')')
end

switch flag
    case 0, data(ind2) = 0; % excluding doublets and above
    case 1,
        data = zeros(size(data)); % including only doublets/triplets/quadruplets
        data(ind) = 1;
end

    function ind = find_doublet_indices(data,dt)
        [a, b] = find(data);
        c = [diff(a); nan;];
        d = find(c>0 & c<dt);
        ind = sub2ind(size(data), a(d), b(d));
    end


end