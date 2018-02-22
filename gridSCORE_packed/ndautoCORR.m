function [acorr,nBins] = ndautoCORR(mat1,mat2,bcut)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This function takes an n-dimensional array and performs an autocorrelation, alternatively two arrays of equal dimensions
%   can be given for a cross-correlation
%   acorr = ndautoCORR(mat1,mat2)
%
%%%%%%%% Inputs
%   mat1 = first array, if only this is provided the function will perform an autocorrelation
%   mat2 = second array, given if we want a cross-correlation, or give a duplicate of mat1 for an autocorrelation
%   bcut = minimum number of bins overlap to include result in autocorrelation
%
%%%%%%%% Outputs
%   acorr = autocorrelation result
%   nBins = the number of bins overlap in autocorrelation, can be used to hide spurious bin values
%
%%%%%%%% Comments
%   12/08/17 adapted from xPearson, changed to be multidimensional, i.e. swap filter2 for imfilter 
%   © Roddy Grieves: rmgrieves@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial variables
if ~exist('mat2','var') || isempty(mat2) || all(isnan(mat2(:)))
    mat2 = mat1;
end
mat1 = double(mat1);
mat2 = double(mat2);

if ~exist('bcut','var') || isempty(bcut) || all(isnan(bcut(:)))
    bcut = 0;
end

filterType = 'full'; % (default = 'full') result type, 'full', 'same' or 'valid'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deal with nans in mat1 or mat2, these should be ignored (i.e. don't contribute to nBins and shouldn't be proliferated by imfilter)
% mat 1
mat1_log = ones(size(mat1));
mat1_log(isnan(mat1)) = 0;
mat1(isnan(mat1)) = 0;
mat1_sqd = mat1.^2;

% mat 2
mat2_log = ones(size(mat2));
mat2_log(isnan(mat2)) = 0;
mat2(isnan(mat2)) = 0;
mat2_sqd = mat2.^2;

%% Prepare vectors for autocorrelation
% autocorrelate two matrices (equivalent to xcorr), basically end up with a matrix as large as both inputs together
mat1_x_mat2 = imfilter(mat1,mat2,filterType);

% calculate bin totals
mat1_sum = imfilter(mat1,mat2_log,filterType); % at each lag, calculate the total of mat1 given the active bins of mat2
mat2_sum = imfilter(mat1_log,mat2,filterType); % at each lag, calculate the total of mat2 given the active bins of mat1
mat1_sqd_sum = imfilter(mat1_sqd,mat2_log,filterType); % at each lag, calculate the total of mat1 squared given the active bins of mat2
mat2_sqd_sum = imfilter(mat1_log,mat2_sqd,filterType); % at each lag, calculate the total of mat2 squared given the active bins of mat1

% at each lag, calculate by how many active bins the two matrices overlap
nBins = imfilter(mat1_log,mat2_log,filterType);
nBins_sqd = nBins.^2; % and find the square of this

%% Perform autocorrelation
mat1_std = (mat1_sqd_sum./nBins - (mat1_sum.^2)./nBins_sqd).^0.5;
mat2_std = (mat2_sqd_sum./nBins - (mat2_sum.^2)./nBins_sqd).^0.5;
covar = mat1_x_mat2 ./ nBins - (mat1_sum.*mat2_sum)./nBins_sqd;
xCoef = covar./(mat1_std.*mat2_std);
acorr = real(xCoef);
acorr(nBins < bcut) = NaN;

if 0
    figure
    subplot(1,2,1)
    mat1(mat1 < 0.1*nanmax(mat1(:))) = NaN;
    vol3d('cdata',mat1,'texture','3D');
    alphamap('rampup');
    % alphamap(0.5 .* alphamap);
    colormap('jet');
    caxis([0 nanmax(mat1(:))]);
    axis on
    view(3)
    daspect([1 1 1])
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis vis3d
    rotate3d on
    title('Matrix in')

    subplot(1,2,2)
    acorr2 = acorr;
    acorr2(acorr2 < 0.2*nanmax(acorr2(:))) = NaN;    
    vol3d('cdata',acorr2,'texture','3D');
    alphamap('rampup');
    colormap('jet');
    caxis([0 nanmax(acorr2(:))]);
    axis on
    view(3)
    daspect([1 1 1])
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis vis3d
    rotate3d on
    title('Autocorr')

end












