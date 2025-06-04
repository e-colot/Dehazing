clear; close all; clc;

image = im2double(imread('./../srcImages/hazy6.png'));

% Ahat = AirlightDirection(image);

% A = airlightAmplitude(image, Ahat);
% disp("End of airlight amplitude estimation");
A = [0.6057; 0.6073; 0.6136];

figure;
imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
title('Estimated Airlight Color A');

[output, transmission] = dehazeHazeLines(image, A);

% denser fog
scaleFactor = 0.5;
denser = scaleFactor*transmission .* output + (1 - scaleFactor*transmission) .* reshape(A, 1, 1, 3);



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
