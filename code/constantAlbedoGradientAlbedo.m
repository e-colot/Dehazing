clear; clc; close all;

src = im2double(imread('../srcImages/hazy1.png'));

A = [0.5 0.5 0.5];

max_itr = 1000;
min_cond = 1e-15;
sensitivity = 1e-4;

for i = 1:max_itr
    vec = zeros(1, 3);
    for k = 1:3
        A_mod = A + 1e-2 * ((1:3)' == k)';
        [t, l] = getParams(src, A_mod);
        C = covariance(t, l)^2;
        vec(k) = C;
    end
    vec = 1./vec;
    vec = vec./norm(vec);
    A = A + sensitivity .* vec;
    [t, l] = getParams(src, A);
    if (covariance(t, l)^2 <= min_cond)
        disp(['Minimizer found']);
        break;
    end
end

[t, l] = getParams(src, A);
tmp = zeros(size(src));
tmp(:,:,1) = (1-t) .* A(1);
tmp(:,:,2) = (1-t) .* A(2);
tmp(:,:,3) = (1-t) .* A(3);
J = (src - tmp) ./ t;

figure;
subplot(121);
imshow(src); title('Hazy Image');
subplot(122);
imshow(J); title('Dehazed Image');


function [transmission, shading] = getParams(img, albedo)
    % In Raanan's paper:
    % transmission is denoted by t
    % shading is denoted by l
    % albedo is denoted by A
    % img is denoted by I

    I_a = (img(:,:,1) .* albedo(1) + img(:,:,2) .* albedo(2) + img(:,:,3) .* albedo(3)) ./ norm(albedo);
    I_r = sqrt((img(:,:,1).^2 + img(:,:,2).^2 + img(:,:,3).^2) - I_a(:,:).^2);

    h = (norm(albedo) - I_a) ./ I_r;

    eta = covariance(I_a, h) / covariance(I_r, h);

    shadingInv = (1 - I_a ./ norm(albedo)) ./ I_r + eta/norm(albedo);
    shading = 1 ./ shadingInv;

    transmission = 1 - (I_a - eta.*I_r) ./ norm(albedo);
end



function [cova] = covariance(v1, v2)

    expextedVal1 = mean(v1, "all");
    expextedVal2 = mean(v2, "all");

    v1 = v1 - expextedVal1;
    v2 = v2 - expextedVal2;

    cova = sum(v1 .* v2, "all") / (size(v1, 1) * size(v1, 2) - 1);

end
