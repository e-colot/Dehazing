clear; close all; clc;

%image = im2double(imread('./../srcImages/hazy1.png'));
image = im2double(imread('./../srcImages/hazy6.png'));

Ahat = AirlightDirection(image);

A = airlightAmplitude(image, Ahat);
%A = Ahat;

figure;
imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
title('Estimated Airlight Color A');

[output, transmission] = dehazeHazeLines(image, A);

% denser fog
scaleFactor = 0.5;
denser = scaleFactor*transmission .* output + (1 - scaleFactor*transmission) .* reshape(A, 1, 1, 3);


% histogram equalization of the output
output = imadjust(output, stretchlim(output), []);

figure;
subplot(221);
imshow(image); 
title('Original Image');
subplot(222);
imshow(output);
title('Dehazed Image');
subplot(223);
imshow(transmission, []);
colormap('jet');
title('Estimated Transmission Map');
subplot(224);
imshow(denser);
title('Denser Fog Simulation');
