clear; close all; clc;

image = im2double(imread('./../srcImages/hazy7.png'));

Ahat = AirlightDirection(image);

A = airlightAmplitude(image, Ahat);
disp("End of airlight amplitude estimation");

% figure;
% imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
% title('Estimated Airlight Color A');

[output, transmission] = dehazeHazeLines(image, A);

% denser fog
scaleFactor = 0.5;
denser = scaleFactor*transmission .* output + (1 - scaleFactor*transmission) .* reshape(A, 1, 1, 3);

stretched = imadjust(output, stretchlim(output), []);

distance = -log(transmission+1e-6);
distance = distance - min(distance(:));
distance = distance / max(distance(:)); % Normalize to [0, 1]

figure;
subplot(121);
imshow(image);
title('Original Image');
subplot(122);
imshow(distance);
colormap('jet');
title('Distance Map');
% figure;
% subplot(221);
% imshow(image); 
% title('Original Image');
% subplot(222);
% imshow(output);
% title('Dehazed Image');
% subplot(223);
% imshow(transmission, []);
% colormap('jet');
% title('Estimated Transmission Map');
% subplot(224);
% imshow(denser);
% title('Denser Fog Simulation');

% figure;
% subplot(131);
% imshow(image);
% title('Input image');
% subplot(132);
% imshow(transmission);
% colormap('jet');
% title('Estimated transmission');
% subplot(133);
% imshow(output);
% title('Dehazed image');


% figure;
% subplot(131);
% imshow(image);
% title('Input image');
% subplot(132);
% imshow(transmission);
% colormap('jet');
% title('Estimated transmission');
% subplot(133);
% imshow(stretched);
% title('Stretched dehazed image');
