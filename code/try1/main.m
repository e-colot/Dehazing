clear; close all; clc;

image = im2double(imread('../../srcImages/hazy3.png'));
output = zeros(size(image));

blockSize = [24, 24];
[rows, cols, ~] = size(image);

for i = 1:blockSize(1):rows
    for j = 1:blockSize(2):cols
        rowEnd = min(i + blockSize(1) - 1, rows);
        colEnd = min(j + blockSize(2) - 1, cols);
        block = image(i:rowEnd, j:colEnd, :);
        dehazedBlock = dehaze(block);
        output(i:rowEnd, j:colEnd, :) = dehazedBlock;
    end
end

figure;
subplot(121);
imshow(image); 
title('Original Image');
subplot(122);
imshow(output);
title('Dehazed Image');
