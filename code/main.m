clear; close all; clc;

%image = im2double(imread('./../srcImages/hazy1.png'));
image = im2double(imread('./../srcImages/hazy3.png'));

Ahat = AirlightDirection(image);

A = airlightAmplitude(image, Ahat);
%A = Ahat;

figure;
imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
title('Estimated Airlight Color A');

[output, transmission] = dehazeHazeLines(image, A);

% histogram equalization of the output
output = imadjust(output, stretchlim(output), []);

figure;
subplot(131);
imshow(image); 
title('Original Image');
subplot(132);
imshow(output);
title('Dehazed Image');
subplot(133);
imshow(transmission, []);
title('Estimated Transmission Map');
