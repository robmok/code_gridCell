function [g,gdata] = gridSCORE(im,method,doPlot) %RM edited 171127 added doPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an autocorrelogram image and calculates various grid metrics
%   [g,gdata] = gridSCORE(im)
%
%%%%%%%% Inputs
%   im = the autocorrelogram image
%
%%%%%%%% Outputs
%   g = grid score
%   gdata = structure containing more detailed info
%         gdata.mid_peak = middle peak location [x,y]
%         gdata.near_peaks = surrounding 6 peak locations [x,y]
%         gdata.near_peaks_d = median distance to surrounding peaks
%         gdata.central_ring = autocorrelogram image cut to central fields (with middle removed)
%         gdata.orientation = the angle/orientation of the grid
%         gdata.g_score = grid score
%
%%%%%%%% Comments
%   21/08/17 created 
%   12/09/17 added Allen and Wills methods
%   13/09/17 added peak method
%   15/10/17 removed buddy peak method and corrected Allen method, added orientation
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('method','var') || isempty(method) || all(isnan(method(:)))
    method = 'allen';
end

% preallocate measures
g = NaN; % gscore
gdata = struct; % structure of analysis details
gdata.mid_peak = NaN;
gdata.near_peaks = NaN;
gdata.near_peaks_d = NaN;
gdata.central_ring = NaN;
gdata.g_score = NaN;
gdata.wavelength = NaN;
gdata.radius = NaN;
gdata.central_mask = NaN;
gdata.mean_inner_angle = NaN;
gdata.orientation = NaN;
gdata.method = method;

im = single(im);
if all(isnan(im(:)))
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate grid score
switch method       
    case {'allen'}
        %% Method used by Perez-Escobar et al. (2016) Visual landmarks sharpen grid cell metric and confer context specificity to neurons of the medial entorhinal cortex
        % find blobs
        imb = im>0.1; %threshold of autocorrelogram
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList'); % get all blobs
        as = [blobs.Area].';
        blobs = blobs(as>10,:); %get blobs more than 10 pixels

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
        cent = [size(im,2)/2,size(im,1)/2];
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % find middle peak and recalculate distance to it instead
        [~,pindx] = min(ds); % find peak closest to centre - this is the origin
        cent = cents(pindx,:);
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % sort blobs according to distance
        [ds,sindx] = sort(ds,'ascend');
        blobs = blobs(sindx,:);
        if length(blobs)==1 %if only 1 blob, exits function
            return
        end
        if length(blobs)>7 %if more than 7, get closest 7
            blobs = blobs(1:7,:);
            ds = ds(1:7,:);
        end       
        
        % calculate grid orientation
        L = createLine(ones(length(blobs)-1,2).*blobs(1).Centroid,cell2mat({blobs(2:end).Centroid}')); %draw line through 6 blob centres, through centre blob
        as = 360 - rad2deg(lineAngle(L));
        grid_ori = min(as); %get smallest angle as grid orientation
        
        % calculate mean distance to closest blobs
        mds = mean(ds);
        dcut = ceil(mds*1.25); %this value determine how much is cut - increasing it to 1.5 --> increased g of the example a bit (from 1.023 to 1.066)
        dcuti = ceil(mds*0.4);

        % cut to the central portion of the autocorrelation - %NOTE this
        % cuts half of each of the 6 peaks out (since takes mean distance
        % of blob centroids to centre, then cuts a circle); will this
        % underestimate?
        rcent = round([blobs(1).Centroid(2),blobs(1).Centroid(1)]);
        imp = padarray(im,[dcut dcut],NaN,'both'); % pad array - sometimes mds is calculated diagonally and is larger than im is wide
        imcent = imp(rcent(1)+dcut-dcut:rcent(1)+dcut+dcut,rcent(2)+dcut-dcut:rcent(2)+dcut+dcut); % take the central part of the padded image
        dmat = zeros(size(imcent));
        dmat(ceil(size(imcent,2)/2),ceil(size(imcent,1)/2)) = 1;
        dmat = bwdist(dmat);
        imcent(dmat>dcut) = NaN;

        % remove the central peak
        imcent(dmat<dcuti) = NaN;

        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(imcent,as(a),'bilinear','crop');
            r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); 
            if ~isnan(r) %RM added 28/02/18 - trapz left/right gridness problem; when >1 blob but <7
                rs(a) = r(1,2);
            else
                rs(a) = nan;
            end
        end
        g = ((rs(2)+rs(4))/2) - ((rs(1)+rs(3)+rs(5))/3);

        % collect data
        gdata.mid_peak = blobs(1).Centroid;
        gdata.near_peaks = cell2mat({blobs(2:end).Centroid}.');
        gdata.near_peaks_d = ds(2:end);
        gdata.central_ring = imcent;

        gdata.g_score = g;
        gdata.wavelength = mds;
        gdata.radius = sqrt(blobs(1).Area)./pi;
        gdata.orientation = grid_ori;

        dmat = zeros(size(im));
        dmat(rcent(1),rcent(2)) = 1;
        dmat = bwdist(dmat);
        msk = ones(size(im)).*0.2;
        msk(dmat<dcut & dmat>dcuti) = 1;
        gdata.central_mask = msk;

        % create figure if required
        if doPlot %1 %RM edited 171127
%             figure
            imc = imagesc(im);
            set(imc,'alphadata',msk);    
            hold on
            plot(gdata.near_peaks(:,1),gdata.near_peaks(:,2),'kx','MarkerSize',10);
            title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));   
            drawLine(L);
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
%             keyboard
        end        
        
    case {'wills'}
        %% Method used by Wills et al. (2012) The abrupt development of adult-like grid cell firing in the medial entorhinal cortex
        % find blobs
        imb = im>0.3;
        blobs = regionprops(imb,'Centroid','Area','PixelIdxList');
        as = [blobs.Area].';
        blobs = blobs(as>3,:);

        % get distance to image centre
        cents = cell2mat({blobs.Centroid}.');
        cent = [size(im,2)/2,size(im,1)/2];
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % find middle peak and recalculate distance to it instead
        [~,pindx] = min(ds); % find peak closest to centre - this is the origin
        cent = cents(pindx,:);
        ds = sqrt(sum((cents-cent).^2,2)); % all distances to centre of image

        % sort blobs according to distance
        [ds,sindx] = sort(ds,'ascend');
        blobs = blobs(sindx,:);
        if length(blobs)==1
            return
        end
        if length(blobs)>7
            blobs = blobs(1:7,:);
            ds = ds(1:7,:);
        end       
        
        % calculate grid orientation
        L = createLine(ones(length(blobs)-1,2).*blobs(1).Centroid,cell2mat({blobs(2:end).Centroid}'));
        as = 360 - rad2deg(lineAngle(L));
        grid_ori = min(as);        
        
        % calculate mean distance to closest blobs
        mds = mean(ds);
        dcut = ceil(mds*1.25);
        dcuti = ceil(mds*0.4);

        % cut to the central portion of the autocorrelation
        rcent = round([blobs(1).Centroid(2),blobs(1).Centroid(1)]);
        imp = padarray(im,[dcut dcut],NaN,'both'); % pad array - sometimes mds is calculated diagonally and is larger than im is wide
        imcent = imp(rcent(1)+dcut-dcut:rcent(1)+dcut+dcut,rcent(2)+dcut-dcut:rcent(2)+dcut+dcut); % take the central part of the padded image
        dmat = zeros(size(imcent));
        dmat(ceil(size(imcent,2)/2),ceil(size(imcent,1)/2)) = 1;
        dmat = bwdist(dmat);
        imcent(dmat>dcut) = NaN;

        % remove the central peak
        imcent(dmat<dcuti) = NaN;

        % rotational correlation
        rs = NaN(5,1);
        as = 30:30:150;
        for a = 1:length(as)
            mrot = imrotate(imcent,as(a),'bilinear','crop');
            r = corrcoef(mrot(:),imcent(:),'rows','pairwise'); 
            if ~isnan(r) %RM added 28/02/18 - trapz left/right gridness problem
                rs(a) = r(1,2);
            else
                rs(a) = nan;
            end
        end
        g = nanmin([rs(2),rs(4)]) - nanmax([rs(1),rs(3),rs(5)]); % big diff to above here: get minimum of the 60/120 deg corr minus highest other

        % collect data
        gdata.mid_peak = blobs(1).Centroid;
        gdata.near_peaks = cell2mat({blobs(2:end).Centroid}.');
        gdata.near_peaks_d = ds(2:end);
        gdata.central_ring = imcent;

        gdata.g_score = g;
        gdata.wavelength = mds;
        gdata.radius = sqrt(blobs(1).Area)./pi;
        gdata.orientation = grid_ori;

        dmat = zeros(size(im));
        dmat(rcent(1),rcent(2)) = 1;
        dmat = bwdist(dmat);
        msk = ones(size(im)).*0.2;
        msk(dmat<dcut & dmat>dcuti) = 1;
        gdata.central_mask = msk;

        % create figure if required
        if doPlot %1 %RM edited 171127
%             figure
            imc = imagesc(im);
            set(imc,'alphadata',msk);    
            hold on
            plot(gdata.near_peaks(:,1),gdata.near_peaks(:,2),'kx','MarkerSize',10);
            title(sprintf('g = %.2f, s = %.2f, r = %.2f, o = %.2f',g,gdata.wavelength,gdata.radius,gdata.orientation));   
            drawLine(L);
            caxis([0 nanmax(imcent(:))])
            daspect([1 1 1]);
            axis off
%             keyboard
        end
end



















        



