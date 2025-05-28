function [J, t] = dehaze(src, A)

    [t, l] = getParams(src, A);
    tmp = zeros(size(src));
    tmp(:,:,1) = (1-t) .* A(1);
    tmp(:,:,2) = (1-t) .* A(2);
    tmp(:,:,3) = (1-t) .* A(3);
    J = (src - tmp) ./ t;

end


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
