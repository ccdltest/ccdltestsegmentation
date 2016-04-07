% BME42-731 / ECE18-795 / CB02-740 Bioimage Informatics
% Spring 2016
% Project Assignment #3
% Due on Apr-04-2016
%
% Ebony Calloway
% Youngsuk Kim
% Alexander Ruesch
% ----------------------------------------------------------------------- %

%% Load data
clear, close all
% general parameter definition
sourcefolder = 'images\';
destinationfolder = 'results\croppedImages\';
nImages = 218;

% Crop background
% Iterate over all images and apply the selected region.
for i=1:nImages
    name = [sourcefolder sprintf('001_a5_002_t%03.0f.tif',i)];
    I = imread(name);
    %     I = I./max(I(:));
    % ask for a sub region only for the first frame.
    if i==1
        [~,pos] = subimg(I);
    end
    % write the image with the name "background001.tif". The choosen format
    % adds leading zeros to the small numbers.
    path = [destinationfolder sprintf('background%03.0f.tif',i)];
    % apply sub region to all images.
    applysubimg(I, pos, path);
end
%% Process data
% B.2.1.
% Calculate mean and standard deviation of background noise.

BG_mu = 0;
BG_sigma = 0;
% Loop over all background images
for i=1:nImages
    name = [destinationfolder sprintf('background%03.0f.tif',i)];
    I = double(imread(name));
    I=I./max(I(:));
    BG_mu = BG_mu + mean(I(:));
    BG_sigma = BG_sigma + std(I(:));
end

BG_mu = BG_mu/nImages;
BG_sigma = BG_sigma/nImages;

%% B.2.2.
% B.2.2 Detection of local maxima and local minima (10 points)
% - Filter each frame with a Gaussian kernel with standard deviation equal
% one third of the Rayleigh radius. The image sequence was collected using
% an objective lens with 100x and a NA of 1.4. The fluorophore used is YFP
% (Yellow Fluorescent Protein). Assume its excitation wavelength at 515 nm.
% Assume a pixel size of 65nm. - Use a 3x3 mask to detect local maxima and
% local minima. - Select one frame; compare detection results using a 3x3
% mask versus a 5x5 mask.

for i=1:nImages
    name = [sourcefolder sprintf('001_a5_002_t%03.0f.tif',i)];
    I = double(imread(name));
    Im{i} = I./max(I(:));
end
Mag = 100;
NA = 1.4;
Wavelength = 515*((10)^-9); %m
PixelSize = 65*((10)^-9); %m

Rr = .61*Wavelength/NA;% radius

EPixelSize = PixelSize/Mag;
%sigma = Rr*EPixelSize/3;
sigma = Rr/3;

k = round(6*sigma) + 1;
%gauss = (1/2*pi*(sigma)^2) * exp(-((gaussx.^2) + (gaussy.^2))/(2*(sigma^2)));
h = fspecial('gaussian', k, sigma);
fImage = imfilter(Im{1},h,'replicate');
subplot(2,1,1)
imshow(Im{1},[])
title('Original Image')
subplot(2,1,2)
imshow(fImage,[])
title('Filtered Image')
pause
[height,width] = size(fImage);
%Local maxima, 3x3
z = 1;
for i=1:height-3
    for j=1:width-3
        mask = fImage(i:i+2,j:j+2);
        if mask(5)== max(mask(:))
            maxima3(z,1) = i+1;
            maxima3(z,2) = j+1;
            z = z + 1;
        end
    end
end
max3 = zeros(size(fImage));
sz =size(maxima3,1);
for i=1:sz
    max3(maxima3(i,1),maxima3(i,2)) = 1;
end

%Local minima, 3x3
z = 1;
for i=1:height-3
    for j=1:width-3
        mask = fImage(i:i+2,j:j+2);
        if mask(5)== min(mask(:))
            minima3(z,1) = i+1;
            minima3(z,2) = j+1;
            z = z + 1;
        end
    end
end
min3 = zeros(size(fImage));
sz =size(minima3,1);
for i=1:sz
    min3(minima3(i,1),minima3(i,2)) = 1;
end
subplot(2,1,1)
imshow(min3,[])
title('Local Minima')
subplot(2,1,2)
imshow(max3,[])
title('Local Maxima')
pause
%Local maxima, 5x5
z = 1;
for i=1:height-5
    for j=1:width-5
        mask = fImage(i:i+4,j:j+4);
        if mask(13)== max(mask(:))
            maxima5(z,1) = i+2;
            maxima5(z,2) = j+2;
            z = z + 1;
        end
    end
end
max5 = zeros(size(fImage));
sz =size(maxima5,1);
for i=1:sz
    max5(maxima5(i,1),maxima5(i,2)) = 1;
end
imshow(max5,[])

%Local minima, 5x5
z = 1;
for i=1:height-5
    for j=1:width-5
        mask = fImage(i:i+4,j:j+4);
        if mask(13)== min(mask(:))
            minima5(z,1) = i+2;
            minima5(z,2) = j+2;
            z = z + 1;
        end
    end
end
min5 = zeros(size(fImage));
sz =size(minima5,1);
for i=1:sz
    min5(minima5(i,1),minima5(i,2)) = 1;
end

subplot(2,2,1)
imshow(min3,[])
title('Local Minima, 3x3')
subplot(2,2,2)
imshow(max3,[])
title('Local Maxima, 3x3')
subplot(2,2,3)
imshow(min5,[])
title('Local Minima, 5x5')
subplot(2,2,4)
imshow(max5,[])
title('Local Maxima, 5x5')
% B.2.3.
%% show data
DT = delaunayTriangulation(minima5(:,2),minima5(:,1));
figure;
imagesc(I); colormap('gray')
hold on;
triplot(DT);
scatter(maxima5(:,2),maxima5(:,1),'r')

%% B.2.4
% Apply statisticel speckle selection from Ponti et al.
Qa = 2.5; % confidence quantile
Mag = 100;
NA = 1.4;
Wavelength = 515*((10)^-9); %m
PixelSize = 65*((10)^-9); %m
Rr = .61*Wavelength/NA;% radius
EPixelSize = PixelSize/Mag;
sigma = Rr/3;
k = round(6*sigma) + 1;
h = fspecial('gaussian', k, sigma);

% BG_mu is taken from B.2.1
% BG_sigma is taken from B.2.1

destinationfolder = 'results\particles\';

% The above shown algorithms are organized in functions for the algorithm
% below. Based on a statistical approach from Ponti et al, we will run the
% particle detection for every frame. The results are stored in speparate
% .mat files on the hard drive.
w=waitbar(0,'Processing time frames...');
for n=1:nImages
    name = [sourcefolder sprintf('001_a5_002_t%03.0f.tif',n)];
    I = double(imread(name));
    I=I./max(I(:));
    I = imfilter(I,h,'replicate');
    
    % Find minima and maxima
    [minima5, maxima5] = localMinMax5(I);
    
    % Perform the triangulation
    DT = delaunayTriangulation(minima5(:,2),minima5(:,1));
    
    for i=1:size(maxima5,1)
        minMaxCorrelation = DT.pointLocation(maxima5(i,2),maxima5(i,1));
        if isnan(minMaxCorrelation)
            maxima5(i,:) = nan;
        else
            corners = DT.ConnectivityList(minMaxCorrelation,:);
            cornerCoord = DT.Points(corners,:);
            
            % calclulate the background noise for every found local maximum.
            I_BG = 0;
            for j=1:3
                I_BG = I_BG + I(cornerCoord(j,2),cornerCoord(j,1));
            end
            I_BG = I_BG / 3;
            deltaI = I(maxima5(i,1),maxima5(i,2)) - I_BG;
            T = abs(deltaI)/(BG_sigma/sqrt(3));
            if T < Qa
                maxima5(i,:) = nan;
            end
        end
    end
    
    particle = maxima5(~any(isnan(maxima5),2),:);
    name = [destinationfolder sprintf('particles%03.0f.mat',n)];
    save(name,'particle');
    waitbar(n/nImages);
end
close(w);

%% Show a slow movie showing the results
disp('To run the video press a random key.')
pause;
figure;
for i=1:nImages
    name = [sourcefolder sprintf('001_a5_002_t%03.0f.tif',i)];
    I = double(imread(name));
    I=I./max(I(:));
    name = [destinationfolder sprintf('particles%03.0f.mat',i)];
    tmp = load(name);
    particle = tmp.particle;
    
    imagesc(1:1392,1:204,I); colormap('gray');
    hold on; 
    scatter(particle(:,2),particle(:,1),'r');
    hold off;
    % Set frame rate here: 
    pause(0.1);
end

 


