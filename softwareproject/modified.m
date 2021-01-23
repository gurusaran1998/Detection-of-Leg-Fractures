clc
clear all;


[filename,pathname] = uigetfile({'*.*';'*.bmp';'*.tif';'*.gif';'*.png';'*.jpg';'*.jpeg'},'Pick an Image File');
img = imread([pathname,filename]);

imgblursigma = 2;   % amount denoiseinput image 

minhoughpeakdistance = 5;   % distance between peaks in hough transform angle detection
 
houghconvolutionlength = 40;        % length of line use detect bone regions

houghconvolutiondilate = 2;             % amount dilate kernel bone detection

breaklinetolerance = 0.25;	% tolerance bone end detection 
 
breakpointdilate = 6;		 % amount dilate detected bone end points

img = (rgb2gray(img)); 		% load image

img = imfilter(img, fspecial('gaussian', 10, imgblursigma), 'symmetric'); % denoise 
 
% edge detection find bone edges in image

% filter out 2 longest lines % feature may need changed if break not in middle of bone
 
boneedges = edge(img, 'canny');

boneedges = bwmorph(boneedges, 'close'); 

edgeregs = regionprops(boneedges, 'area', 'pixelidxlist'); 

arealist = sort(vertcat(edgeregs.area), 'descend'); 

edgeregs(~ismember(vertcat(edgeregs.area), arealist(1:2))) = []; 

edgeimg = zeros(size(img, 1), size(img,2)); 

edgeimg(vertcat(edgeregs.pixelidxlist)) = 1;  

% hough transform on edge image find angles @ bone pieces 
% found 
% use max value of hough transform vs angle find angles @ lines 
% oriented.  if there more 1 major angle contribution there 
% 2 peaks detected 1 peak if there 1 major % angle contribution (ie peaks here = number of located bones = number of % breaks + 1) 

[h,t,r] = hough(edgeimg,'rhoresolution',1,'theta',-90:2:89.5); 

maxhough = max(h, [], 1); 

houghthresh = (max(maxhough) - min(maxhough))/2 + min(maxhough); 

[~, houghpeaks] = findpeaks(maxhough,'minpeakheight',houghthresh, 'minpeakdistance', minhoughpeakdistance);
  
% plot hough detection results 

figure(1) 

plot(t, maxhough); 

hold on 

plot([min(t) max(t)], [houghthresh, houghthresh], 'r'); 

plot(t(houghpeaks), maxhough(houghpeaks), 'rx', 'markersize', 12, 'linewidth', 2); 

hold off 

xlabel('theta value'); 

ylabel('max hough transform'); 

legend({'max hough transform', 'hough peak threshold', 'detected peak'});

% locate site of break 

if numel(houghpeaks) > 1;  

    breakstack = zeros(size(img, 1), size(img, 2), numel(houghpeaks)); 
    
    % convolute edge image line of detected angle hough transform 
    for m = 1:numel(houghpeaks);

        bonekernel = strel('line', houghconvolutionlength, t(houghpeaks(m)));
        
        kern = double(bwmorph(bonekernel.getnhood(), 'dilate', houghconvolutiondilate));
        
        breakstack(:,:,m) = imfilter(edgeimg, kern).*edgeimg;
        
        
        
    end

% take difference between convolution images.  crosses 0
  
% (within tolerance) should break is.  have filter out   
  
% regions elsewhere bone ends.
    
brimg = abs(diff(breakstack, 1, 3)) < breaklinetolerance*max(breakstack(:)) & edgeimg > 0;  
   
[bpy, bpx] = find(abs(diff(breakstack, 1, 3)) < breaklinetolerance*max(breakstack(:)) & edgeimg > 0); 
    
brimg = bwmorph(brimg, 'dilate', breakpointdilate);   

brreg = regionprops(brimg, 'area', 'majoraxislength', 'minoraxislength', ... 
        'orientation', 'centroid');   
  
brreg(vertcat(brreg.area) ~= max(vertcat(brreg.area))) = [];  
    
% calculate bounding ellipse  
  
brreg.ellipsecoords = zeros(100, 2);  
   
t = linspace(0, 2*pi, 100);  
   
brreg.ellipsecoords(:,1) = brreg.centroid(1) + brreg.majoraxislength/2*cos(t - brreg.orientation);     brreg.ellipsecoords(:,2) = brreg.centroid(2) + brreg.minoraxislength/2*sin(t - brreg.orientation);  

else     
brreg = []; 
 
end

% draw ellipse around break location
 
figure(2)
 
imshow(img) 

hold on 

colormap('gray') 

if ~isempty(brreg)  
   
    plot(brreg.ellipsecoords(:,1), brreg.ellipsecoords(:,2), 'r');
    
    msgbox('Bone Fracture Detected')
    
else
    
    msgbox('Normal, No Fracture Bone')
 
end

hold off
