clear; close all; clc;

image = im2double(imread('./../srcImages/hazy4.jpg'));

Ahat = AirlightDirection(image);

A = airlightAmplitude(image, Ahat);
disp("End of airlight amplitude estimation");

% figure;
% imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
% title('Estimated Airlight Color A');
A = Ahat;
[output, transmission] = dehazeHazeLines(image, A);

% denser fog
scaleFactor = 0.5;
denser = scaleFactor*transmission .* output + (1 - scaleFactor*transmission) .* reshape(A, 1, 1, 3);

stretched = imadjust(output, stretchlim(output), []);

distance = -log(transmission+1e-6);
distance = distance - min(distance(:));
distance = distance / max(distance(:)); % Normalize to [0, 1]

figure;
subplot(1, 3, 1);
imshow(image);
title('Original Hazy Image');
subplot(1, 3, 2);
imshow(output);
title('Dehazed Image');
subplot(1, 3, 3);
imshow(transmission);
colormap("jet");
title('Estimated Transmission Map');
