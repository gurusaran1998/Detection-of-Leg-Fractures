clc
clear all;


[filename,pathname] = uigetfile({'*.*';'*.bmp';'*.tif';'*.gif';'*.png';'*.jpg';'*.jpeg'},'Pick an Image File');
img = imread([pathname,filename]);

figure(1)

imshow(img); % distance between peaks in hough transform angle detection

title('Input Image');

ImgBlurSigma = 2;              % amount denoiseinput image 

MinHoughPeakDistance = 5;  % distance between peaks in hough transform angle detection

HoughConvolutionLength = 40; % Length of line to use to  detect bone regions

HoughConvolutionDilate = 2;  % amount dilate kernel bone detection

BreakLineTolerance = 0.25;   % tolerance bone end detection

breakPointDilate = 6;        % amount dilate detected bone end points

img_gray = (rgb2gray(img)); 

img_filtered = imfilter(img_gray, fspecial('gaussian', 10, ImgBlurSigma), 'symmetric'); % denoise 

figure(2)

imshow(img_filtered);

title('Filtered image');

% edge detection find bone edges in image

% filter out 2 longest lines % feature may need changed if break not in middle of bone

boneEdges = edge(img_filtered, 'canny');

boneEdges1 = bwmorph(boneEdges, 'close');

figure(3)

imshow(boneEdges1);

title('Edge Detection');

edgeRegs = regionprops(boneEdges1, 'Area', 'PixelIdxList');

AreaList = sort(vertcat(edgeRegs.Area), 'descend');

edgeRegs(~ismember(vertcat(edgeRegs.Area), AreaList(1:2))) = [];

edgeImg = zeros(size(img_filtered, 1), size(img_filtered,2));

edgeImg(vertcat(edgeRegs.PixelIdxList)) = 1;
 
% hough transform on edge image find angles @ bone pieces 
% found 
% use max value of hough transform vs angle find angles @ lines 
% oriented.  if there more 1 major angle contribution there 
% 2 peaks detected 1 peak if there 1 major % angle contribution (ie peaks here = number of located bones = number of % breaks + 1) 

[H,T,R] = hough(edgeImg,'RhoResolution',1,'Theta',-90:2:89.5);

maxHough = max(H, [], 1);

HoughThresh = (max(maxHough) - min(maxHough))/2 + min(maxHough);

[~, HoughPeaks] = findpeaks(maxHough,'MINPEAKHEIGHT',HoughThresh, 'MinPeakDistance', MinHoughPeakDistance);


% plot hough detection results 

figure(1) 

plot(T, maxHough); 

hold on 

plot([min(T) max(T)], [HoughThresh, HoughThresh], 'R'); 

plot(T(HoughPeaks), maxHough(HoughPeaks), 'rx', 'markersize', 12, 'linewidth', 2); 

hold off 

xlabel('theta value'); 

ylabel('max hough transform'); 

legend({'max hough transform', 'hough peak threshold', 'detected peak'});



% locate site of break 
if numel(HoughPeaks) > 1;
    
    BreakStack = zeros(size(img_filtered, 1), size(img_filtered, 2), numel(HoughPeaks));
    
    % Convolute edge image with line of detected angle from hough transform
    
    for m = 1:numel(HoughPeaks);

        boneKernel = strel('line', HoughConvolutionLength, T(HoughPeaks(m)));
        
        kern = double(bwmorph(boneKernel.getnhood(), 'dilate', HoughConvolutionDilate));
        
        BreakStack(:,:,m) = imfilter(edgeImg, kern).*edgeImg;
        
        
        
    end
    
% take difference between convolution images.  crosses 0
  
% (within tolerance) should break is.  have filter out   
  
% regions elsewhere bone ends.
        
    brImg = abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0;
    
    [BpY, BpX] = find(abs(diff(BreakStack, 1, 3)) < BreakLineTolerance*max(BreakStack(:)) & edgeImg > 0);
    
    brImg = bwmorph(brImg, 'dilate', breakPointDilate);
    
    figure(5);
    
    imshow(brImg);
    
    brReg = regionprops(brImg, 'Area', 'MajorAxisLength', 'MinorAxisLength', ...
        'Orientation', 'Centroid');
    
    brReg(vertcat(brReg.Area) ~= max(vertcat(brReg.Area))) = [];

    % Calculate bounding
    
    brReg.EllipseCoords = zeros(100, 2);
    
    t = linspace(0, 2*pi, 100);
    
    brReg.EllipseCoords(:,1) = brReg.Centroid(1) + brReg.MajorAxisLength/2*cos(t - brReg.Orientation);
    
    brReg.EllipseCoords(:,2) = brReg.Centroid(2) + brReg.MinorAxisLength/2*sin(t - brReg.Orientation);

else
    
    brReg = [];

end

% draw ellipse around break location

figure(6)

imshow(img)

title('Detected Output');

hold on

colormap('gray')

if ~isempty(brReg)
    
    plot(brReg.EllipseCoords(:,1), brReg.EllipseCoords(:,2), 'r');
    
    msgbox('Bone Fracture Detected')

else
    
    msgbox('Normal, No Fracture Bone')

end

hold off
    