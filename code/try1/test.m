clear; close all; clc;



NS = 64 ;
N = NS * NS ;

A = [0.8 0.8 0.9]' ;

R = [0.8 0.4 0.8]' ;

Rt = R - A * (R' * A) / (norm(A)^2);

true_eta = A' * R / (norm(A) * norm(Rt)) ;

%
% Generating synthetic transmission & shading functions
%

i=0;
for y=1:NS
    for x=1:NS
    i=i+1 ;
    
    xx = x / NS - 0.5 ;
    yy = y / NS - 0.5 ;

    rr = norm([xx yy]) ;
    
    t(i) = 0.3 + 0.2 * sin(rr * 60) + 0.1 * rand ;
    l(i) = 1 + 0.3 * sin((xx+yy) * 80) + 0.1 * rand ;
    
    end
end

t=t';
l=l';
    

%
% Generating synthetic input image
%

i=0;
for y=1:NS
    for x=1:NS
    i=i+1 ;

    I(i,:) = t(i) .* l(i) .* R + (1-t(i)) .* A ;
    
    im(x,y,1) = I(i,1) ; im(x,y,2) = I(i,2) ; im(x,y,3) = I(i,3) ;
    tim(x,y,1) = l(i) .* R(1) ; tim(x,y,2) = l(i) .* R(2) ; tim(x,y,3) = l(i) .* R(3) ; 
    end
end

% estimating eta, trans. (and shading)
I = reshape(I, NS, NS, 3);

disp(['Size of I: ', num2str(size(I))]);

[t, l] = getParams(I, A);
tmp = zeros(size(I));
tmp(:,:,1) = (1-t) .* A(1);
tmp(:,:,2) = (1-t) .* A(2);
tmp(:,:,3) = (1-t) .* A(3);
J = (I - tmp) ./ t;

figure;
subplot(121);
imshow(I); title('Hazy Image');
subplot(122);
imshow(J); title('Dehazed Image');





function [transmission, shading] = getParams(img, albedo)
    % In Raanan's paper:
    % transmission is denoted by t
    % shading is denoted by l
    % albedo is denoted by A
    % img is denoted by I

    dimensions = size(img);
    img = reshape(img, [], 3); % Reshape the image to a 2D array

    I_a = img * albedo / norm(albedo);
    I_r = sqrt(sum(img.^2, 2) - I_a.^2);

    h = (norm(albedo) - I_a) ./ I_r;

    eta = covariance(I_a, h) / covariance(I_r, h);

    transmission = 1 - (I_a - eta*I_r) / norm(albedo);
    shading = I_r ./ transmission;

    transmission = reshape(transmission, dimensions(1), dimensions(2));
    shading = reshape(shading, dimensions(1), dimensions(2));
end



function [cova] = covariance(v1, v2)

    expextedVal1 = mean(v1, "all");
    expextedVal2 = mean(v2, "all");

    v1 = v1 - expextedVal1;
    v2 = v2 - expextedVal2;

    cova = sum(v1 .* v2, "all") / (size(v1, 1) * size(v1, 2) - 1);

end