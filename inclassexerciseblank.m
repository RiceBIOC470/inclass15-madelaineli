%% GB comments
Step1 100
Step2 100
Step3 100 
Step4 100
Step5 100
Step6 100 
Step7 100 
Step8 100
Overall 100


%% step 1: write a few lines of code or use FIJI to separately save the
% nuclear channel of the image Colony1.tif for segmentation in Ilastik

img = ('48hColony1_DAPI.tif');
reader = bfGetReader(img);

iplane = reader.getIndex(0,0,0)+1;
plane = bfGetPlane(reader,iplane);
figure(1)
imshow(plane,[]);

%% step 2: train a classifier on the nuclei
% try to get the get nuclei completely but separe them where you can
% save as both simple segmentation and probabilities

% saved as Probability.h5 and Segmentation.h5
%% step 3: use h5read to read your Ilastik simple segmentation
% and display the binary masks produced by Ilastik 

% (datasetname = '/exported_data')
% Ilastik has the image transposed relative to matlab
% values are integers corresponding to segmentation classes you defined,
% figure out which value corresponds to nuclei

immask = h5read('Segmentation.h5', '/exported_data');
immask = squeeze(immask);
mask_seg = immask >= 2;
for ii = 1:size(mask_seg,3)
    mask_seg(:,:,ii) = mask_seg(:,:,ii)';
end
figure(2)
imshow(mask_seg);
%% step 3.1: show segmentation as overlay on raw data
img = imread(img);
figure(3)
img_fuse = imfuse(img,mask_seg,'blend','Scaling','joint');
imshow(img_fuse)
%% step 4: visualize the connected components using label2rgb
% probably a lot of nuclei will be connected into large objects
img_nuc = label2rgb(mask_seg);
figure(4)
imshow(img_nuc)
%% step 5: use h5read to read your Ilastik probabilities and visualize
im_prob = h5read('Probability.h5', '/exported_data');
im_prob = squeeze(im_prob(1,:,:));
figure(5)
imshow(im_prob)
% it will have a channel for each segmentation class you defined

%% step 6: threshold probabilities to separate nuclei better
avr = sum(sum(im_prob))/length(im_prob);
im_prob = im_prob >= (avr-1310);
figure(6)
imshow(im_prob)
%% step 7: watershed to fill in the original segmentation (~hysteresis threshold)
CC = bwconncomp(mask_seg);
stats = regionprops(CC,'Area');
area = [stats.Area];
fusedCandidates = area > mean(area) + std(area);
sublist = CC.PixelIdxList(fusedCandidates);
sublist = cat(1,sublist{:});
fusedMask = false(size(mask_seg));
fusedMask(sublist) = 1;
D = bwdist(~fusedMask);
D = -D;
D(~fusedMask) = -Inf;
L = watershed(D);
rgb = label2rgb(L,'jet',[0.5 0.5 0.5]);
newNuclearMask = L>1 | (mask_seg - fusedMask);
figure(7)
imshow(newNuclearMask)
%% step 8: perform hysteresis thresholding in Ilastik and compare the results
% explain the differences
hyst=h5read('48hColony1_DAPI_Object Predictions.h5', '/exported_data');
hyst=squeeze(hyst);
hyst_rgb=label2rgb(hyst,'jet',[0.5 0.5 0.5]);
figure (8)
imshow(hyst_rgb);
% In hysteresis, cells are more separate than in watershed. Hysteresis
% image has more open spaces in between the segments, whereas in
% watershed, the segments are not as separated.
%% step 9: clean up the results more if you have time 
% using bwmorph, imopen, imclose etc
BW2 = bwmorph(mask_seg,'open');
figure(9)
imshow(BW2);
