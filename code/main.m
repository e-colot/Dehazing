clear; close all; clc;

image = im2double(imread('./../srcImages/hazy1.png'));

Ahat = AirlightDirection(image);

A = airlightAmplitude(image, Ahat);


[output, transmission] = dehazeHazeLines(image, A);

% denser fog
scaleFactor = 0.5;
denser = scaleFactor*transmission .* output + (1 - scaleFactor*transmission) .* reshape(A, 1, 1, 3);

stretched = imadjust(output, stretchlim(output), []);

distance = -log(transmission+1e-6);
distance = distance - min(distance(:));
distance = distance / max(distance(:)); % Normalize to [0, 1]

figure;
subplot(2, 2, 1);
imshow(image);
title('Original Hazy Image');
subplot(2, 2, 2);
imshow(output);
title('Dehazed Image');
subplot(2, 2, 3);
imshow(transmission);
colormap("jet");
title('Estimated Transmission Map');
subplot(2, 2, 4);
imshow(distance);
colormap("jet");
title('Distance Map');
