function diag_autocorr(ensembles)
% Plot autocorrelation per walker and parameter (depends on economic toolbox!)
%
%CUSTOM FUNCTIONS
% acf.m

[npar, nwalk, lchain] = size(ensembles);

% control the number of displayed walers and parameters
nwalk_  = min(nwalk,10);
npar_   = min(npar,5);
if nwalk_ < nwalk
    warning('MATLAB:diagnosisPlot',...
        ['Number of walkers displayed on the autocorrelation plot is cropped to ' num2str(nwalk_)])
end
if npar_ < npar
    warning('MATLAB:diagnosisPlot',...
        ['Number of parameters displayed on the autocorrelation plot is cropped to ' num2str(npar_)])
end

figure('Position', [200, 100, 150*nwalk_, 200*npar_])
for jj = 1:nwalk_
    for ii = 1:npar_
        subplot(npar_, nwalk_, jj+(ii-1)*nwalk_)
%         autocorr(ensembles(ii,jj,:),min(50,lchain-1));
        y = ensembles(ii,jj,:);
        acf(y(:),min(50,lchain-1));
        ylabel(gca,'');
        title(gca,'');
        xlabel(['\theta_{', num2str(ii),'}, walker_{', num2str(jj),'}'])
    end
end


axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Autocorrelation plots','HorizontalAlignment',...
    'center','VerticalAlignment', 'top')

end
