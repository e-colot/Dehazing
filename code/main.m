clear; close all; clc;

image = im2double(imread('./../srcImages/hazy1.png'));

Ahat = AirlightDirection(image);
A = airlightAmplitude(image, Ahat);

figure;
imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
title('Estimated Airlight Color A');

[output, ~] = dehaze(image, A);

figure;
subplot(121);
imshow(image); 
title('Original Image');
subplot(122);
imshow(output);
title('Dehazed Image');
