function diag_triangle(ensembles)
% Plot triangle diagnostic plots
% marginal histograms + pairwise 2D marginals + correlation estimates
% + credible intervals (1D and 2D)
%
%inspired by 
% (1) Hartig et al. (2014). Technical Note: Approximate Bayesian parameterization of a process-based tropical forest mode
% doi:10.5194/bg-11-1261-2014 Fig.3
% (2) triangle.py https://github.com/dfm/triangle.py
%
% TODO: circles instead of square scatter plot marker
%
% NOTE:
%  Due to transparent patches the axes can partially disappear on plots..
%  It is a known bug in Matlab, for details, workaround see:
%  http://stackoverflow.com/questions/9775020/black-lines-missing-in-the-box-holding-the-axes-of-a-matlab-plot
%
%CUSTOM FUNCTIONS
% subtightplot.m
% plot_squares.m
% inpoly.m
% histcn.m
%  Slight modifcation is needed to get the same output format as of
%   _hist_ and _hist3_:
%  the max element is the end of the last bin and counted as being outside of it,
%  thus the _count_ output is larger than _mid_ or _edges_ outputs;
%  putting the max element into the last bin + deleting the last element(s)
%  of _count_ yields tp the desired output format
%

% level of credible interval!
credi_mass = 0.90;

[npar, nwalk, lchain] = size(ensembles);
% squeeze
squeezed_ens = ensembles(:,:)';

npar_   = min(npar,5);
if npar_ < npar
    warning('MATLAB:diagnosisPlot',...
        ['Number of parameters displayed on the ''complex'' histogram plot is cropped to ' num2str(npar_)])
end

figure('Position', [200, 100, 200*npar_, 200*npar_])

% get the histograms 'coodinates' and axis' limits
nbins   = 20;
h       = nan(npar_,nbins);
xh      = nan(npar_,nbins);
xmax    = nan(npar_,1);
xmin    = nan(npar_,1);
%loop over the parameters,
for ii = 1:npar_
    [count, ~, mid] = histcn(squeezed_ens(:,ii),nbins);
    xmax(ii) = max(mid{1});
    xmin(ii) = min(mid{1});
    
    count(end-1)= count(end-1) + count(end);
    count(end)  = [];
    h(ii,:)     = count;
    xh(ii,:)    = mid{1};
end

%PLOT
%loop over the rows
for ii = 1:npar_
    %loop over the columns
    for jj = 1:npar_%ii
        p = sub2ind([npar_, npar_], jj, ii);
        subtightplot(npar_,npar_,p, [0.01, 0.01], 0.1,0.1)

        if ii == jj
            % HISTOGRAM
%             face_color = [50,150,255]/255;
            face_color = [0.4444    0.4913    0.5694];
            bar(xh(ii,:), h(ii,:), 'FaceColor', face_color, 'EdgeColor', [1,1,1])
            set(gca,'YTickLabel','')
            xlim([xmin(jj),xmax(jj)])
            
            % CREDIBLE INTERVAL
            ci = credi_interval(squeezed_ens(:,ii), credi_mass);
            ypos = 0.02*max(h(ii,:));

            hold on
            plot([ci(1), ci(2)], [ypos, ypos],'Color','red','LineWidth',2)
            
        elseif jj < ii
            % pairwise scatter
            % plot(squeezed_ens(:,ii), squeezed_ens(:,jj), '-x')
            
            % 2D HISTOGRAM COUNTING
%             [N0, C0] = hist3(squeezed_ens(:,[jj,ii]),[20,20]);
            [count, ~, mid] = histcn(squeezed_ens(:,[jj,ii]),20,20);
            % the max element is the end of the last bin and counted as being outside of it
            % putting it to the last bin; + get consistent output
            count(end-1,:)= count(end-1,:) + count(end,:);
            count(end,:)  = [];
            count(:,end-1)= count(:,end-1) + count(:,end);
            count(:,end)  = [];
            
            N  = count;
            C  = mid;
            
            % CONTOUR (with default/automatic number of levels)
            CM  = contour(C{1}, C{2}, N');
            cm  = flipud(colormap('bone'));
            ncm = size(cm,1);
            % get the contour lines' coordinates
            cstruct = splitcontours(CM);
            npatch  = size(cstruct,2);
            % get color for areas between contours
            patch_cm = interp1(1:ncm, cm, (1:npatch)/npatch*ncm);

            % FILL AREAS BETWEEN CONTOURS (contourf is harder to set up correctly for me..)
            for kk = 1:npatch
                % color spec
                C = ones(size(cstruct(kk).x.')); 
                patch(cstruct(kk).x.',cstruct(kk).y.',C, 'FaceColor', patch_cm(kk,:))
            end
            
            % ADD POINTS
            % most outer contour for _contour_
            oCM     = [cstruct(1).x; cstruct(1).y];
            % check if the simulation points are inside the outer contour line
            in      = inpoly(squeezed_ens(:,[jj,ii]),oCM');
            hold on
            % plot simulation points outside the outer contour line
            xx      = squeezed_ens(~in,jj);
            yy      = squeezed_ens(~in,ii);
%             plot(xx,yy,'.',...
%                 'MarkerEdgeColor',[0.7, 0.7, 0.7],...
%                 'MarkerSize',1)
            plot_squares(xx,yy,15)
            
            % CREDIBLE CONTOUR

            CI = credi_contour(squeezed_ens(:,[jj,ii]), credi_mass);
            plot(CI(:,1),CI(:,2),'Color','red','LineWidth',2)
            
            % fit straigth line to data
            %             pp   = polyfit(squeezed_ens(:,ii), squeezed_ens(:,jj),1);
            %
            %             xmin    = min(squeezed_ens(:,ii));
            %             xmax    = max(squeezed_ens(:,ii));
            %             dx      = xmax-xmin;
            %             xx      = linspace(xmin+0.1*dx, xmax-0.1*dx,20);
            %             yy      = polyval(pp,xx);
            %             hold on
            %             plot(xx,yy,'Color','black', 'Linewidth',1.5)
            xlim([xmin(jj),xmax(jj)])
            ylim([xmin(ii),xmax(ii)])
            
        else
            % correlation value and coloring
            axis off
            %             rho = corr(squeezed_ens(:,ii), squeezed_ens(:,jj), 'type', 'spearman');
            rho = correlation(squeezed_ens(:,ii), squeezed_ens(:,jj),'spearman');
            
            % color the subplot 'proportional' to _rho_
            %             cmap  = colormap('jet');
            %             c_idx = min(max(round(0.5*size(cmap,1)*(1 + rho)),1), size(cmap,1));
            %                 cmap  = colormap('gray');
            %                 cmap  = flipud(cmap);
            %                 cmap(:,1) = 1;
            %                 c_idx = min(max(round(size(cmap,1)*abs(rho)),1), size(cmap,1));
            %                 facecolor = cmap(c_idx,:);
            %                 patch([0,0,1,1], [0,1,1,0],1, 'FaceColor', facecolor, 'FaceAlpha', 0.5)
            
            % fontsize is proportional to _abs(rho)_
            fontsize = interp1([0,1],[8,16],abs(rho));
            text(0.5,0.5, num2str(round(rho*100)/100), 'FontSize', fontsize,...
                'HorizontalAlignment','center')
        end
        % set axis ticks, x
        if ii == npar_
            %leave the XTick labels
        else
            set(gca,'XTickLabel','')
        end
        % set axis ticks, y
        if jj == 1 && ii > 1
            %leave the YTick labels
        else
            set(gca,'YTickLabel','')
        end
        
        % add axis label
        if jj == 1
            ylabel(['\theta_', num2str(ii)])
        end
        if ii == npar_
            xlabel(['\theta_', num2str(jj)])
        end
        
        % to tackle the disappearing axes (not effective...)
        set(gca,'Layer','top')
    end
end

% add title
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
    'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
if npar_ == 1
    triangle_title = 'Histogram';
elseif npar_ == 2
    triangle_title = '1D(diagonal) marginals and 2D distribution(contour)';
else
    triangle_title = '1D(diagonal) and 2D marginals(contour)';
end
text(0.5, 1,{['\bf', triangle_title, ' with Spearman correlation'],...
    ['and credible interval (hd, level=',num2str(credi_mass),'), p(\theta|D)']},...
    'HorizontalAlignment','center',...
    'VerticalAlignment', 'top')

end