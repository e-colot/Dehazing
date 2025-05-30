clear; close all; clc;

%image = im2double(imread('./../srcImages/hazy1.png'));
image = im2double(imread('./../srcImages/hazy4.jpg'));

Ahat = AirlightDirection(image);

%A = airlightAmplitude(image, Ahat);
A = 0.9*Ahat;
%A = [0.84; 0.79; 0.84]; % Manually set airlight for testing

figure;
imshow(reshape(A, 1, 1, 3), InitialMagnification=10000);
title('Estimated Airlight Color A');

[output, transmission] = dehazeHazeLines(image, A);


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
