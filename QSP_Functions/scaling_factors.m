function c = scaling_factors(scaling,baselineparameters,n_G,n_p,n_V)
scalings.G = scaling(1:length(n_G));
scalings.p = scaling(length(n_G)+1:length(n_p)+length(n_G));
scalings.V = scaling(length(n_G)+length(n_p)+1:end);

for iF = 1:length(n_G)
    aF = n_G{iF};
    c.G.(aF) = baselineparameters.G.(aF) * scalings.G(iF);
end

for iF = 1:length(n_p)
    aF = n_p{iF};
    c.p.(aF) = baselineparameters.p.(aF) * scalings.p(iF);
end

for iF = 1:length(n_V)
    aF = n_V{iF};
    c.V.(aF) = scalings.V(iF) + baselineparameters.V.(aF);
end
