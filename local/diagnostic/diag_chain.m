function diag_chain(ensembles)
% Plot 1D chains (trace)

[npar, nwalk, lchain] = size(ensembles);
% squeeze
squeezed_ens = ensembles(:,:)';

npar_   = min(npar,5);
if npar_ < npar
    warning('MATLAB:diagnosisPlot',...
        ['Number of parameters displayed on the chain plot is cropped to ' num2str(npar_)])
end

figure('Position', [200, 100, 600, 150*npar_])
for ii = 1:npar_
    subplot(npar_, 1, ii)
    plot(squeezed_ens(:,ii))
    xlabel(['\theta_', num2str(ii)])
    xlim([0, size(squeezed_ens,1)])
    if ii == 1
        title('\bf 1D chain (trace) per parameters')
    end
end
end