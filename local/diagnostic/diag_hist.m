function diag_hist(ensembles)
% 1D marginal histogram of chains per parameters
%
%
%CUSTOM FUNCTIONS
% histcn.m

% level of credible interval!
credi_mass = 0.9;

[npar, nwalk, lchain] = size(ensembles);
% squeeze
squeezed_ens = ensembles(:,:)';

npar_   = min(npar,5);
if npar_ < npar
    warning('MATLAB:diagnosisPlot',...
        ['Number of parameters displayed on the histogram plot is cropped to ' num2str(npar_)])
end

figure('Position', [200, 100, 500, 200*npar_])
%loop over the rows
for ii = 1:npar_
    sample_ii = squeezed_ens(:,ii);
    % highes density credible interval
    ci = credi_interval(sample_ii,credi_mass,'hd');
    % histogram
    [count, ~, mid] = histcn(sample_ii,30);
    count(end-1)    = count(end-1) + count(end);
    count(end)      = [];
    mid             = mid{1};
    
    % color differently the mass within the credible interval
    tail_idx = mid < ci(1) | mid > ci(2);
    
    subplot(npar_, 1, ii)
    bar(mid(~tail_idx),count(~tail_idx), 'FaceColor', [200,0,0]/255, 'EdgeColor', 'none')
    hold on
    bar(mid(tail_idx),count(tail_idx), 'FaceColor', [50,150,255]/255, 'EdgeColor', 'none')
    set(gca,'YTickLabel','')
    xlabel(['\theta_', num2str(ii)])
    
    % display median, mean and credible interval on the plot
    ylim = get(gca, 'YLim');
    ypos = ylim(2) - diff(ylim)*0.15;
    xpos = mean(get(gca, 'XLim'));
    str = ['median = ',num2str(roundsd(median(sample_ii),3)),...
        '; mean = ',num2str(roundsd(mean(sample_ii),3)),...
        '; CI_{',num2str(credi_mass),'} = [',num2str(roundsd(ci(1),3)),...
        ', ',num2str(roundsd(ci(2),3)),']'];
    text(xpos,ypos, str, 'HorizontalAlignment', 'center')
end

% add title
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
str = {'\bf Marginal, posterior distibution of the parameters', ['with credible interval(hd, level=',num2str(credi_mass),'),  p(\theta|D)']};
text(0.5, 1, str, 'HorizontalAlignment', 'center','VerticalAlignment', 'top')
end