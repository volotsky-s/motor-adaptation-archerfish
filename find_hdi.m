function hdi_lim = find_hdi(sample_vec,cred_mass)

%   sampleVec
%     vector - representative values from a probability distribution
%   credMass
%     scalar between 0 and 1 - indicating the mass within the credible
%     interval that is to be estimated
%
%   hdiLim - vector - limits of HDI

sorted_vec = sort(sample_vec);
inds_hdi = ceil(cred_mass * length(sorted_vec));
n_inds_hdi = length(sorted_vec) - inds_hdi; 

int_width = zeros(n_inds_hdi,1); 
for ind = 1:n_inds_hdi 
    int_width(ind) = sorted_vec(ind + inds_hdi) - sorted_vec(ind);
end

[~,idx_min] = min(int_width);
hdi_low_lim = sorted_vec(idx_min);
hdi_high_lim = sorted_vec(idx_min + inds_hdi);
hdi_lim = [hdi_low_lim, hdi_high_lim];

end
