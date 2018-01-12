function mapping = calc_projection_from_library(endmembers, reduced_dim)
if 1 % Make the same size for all endmembers and the mean lie at the center
    num_spectra_per_endm = zeros(1,length(endmembers));
    for i = 1:length(endmembers)
        num_spectra_per_endm(i) = size(endmembers{i},1);
    end
    num_spectra = round(mean(num_spectra_per_endm));
    spectra_train = [];
    for i = 1:length(endmembers)
        inds = linspace(1, num_spectra_per_endm(i), num_spectra);
        inds = round(inds);
        spectra_train = cat(1, spectra_train, endmembers{i}(inds,:));
    end
else % The size difference is usually not big for spectra libraries.
    spectra_train = cell2mat(endmembers');
end
[~,mapping] = pca(spectra_train, reduced_dim);


end

