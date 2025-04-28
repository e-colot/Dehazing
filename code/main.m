clear; close all; clc;

image = im2double(imread('../../srcImages/hazy1.png'));
output = dehaze(image);

figure;
subplot(121);
imshow(image); 
title('Original Image');
subplot(122);
imshow(output);
title('Dehazed Image');
