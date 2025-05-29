function [J, t_est] = dehazeHazeLines(I, A)
% Dehazing using Haze-Lines (Berman et al., CVPR 2016)
% I: input hazy image (HxWx3)
% A: airlight (1x3)
% J: dehazed image
% t: estimated transmission map

% Step 1: Subtract airlight
    [h, w, ~] = size(I);
    IA = double(I) - reshape(A, 1, 1, 3);

% Step 2: Convert to spherical coordinates
    [R, Phi, Theta] = rgb2spherical(IA);

% Step 3: define haze lines by clustering pixels [Phi, Theta]
    clusterCnt = 1000;
    X = [Phi(:), Theta(:)];
    [idx, ~] = kmeans(X, clusterCnt, 'MaxIter', 200);

% Prepare output
t_est = zeros(size(R));

% Step 4: For each haze-line
    for H = 1:clusterCnt
        mask = reshape(idx == H, size(R));
        rInH = R(mask);
        % Step 5: Estimate max radius in cluster H
        r_max = max(rInH);
        % Step 6+7: For each pixel in cluster, estimate transmission
        t_est(mask) = rInH / r_max;
    end

% Step 8: regularization

    % Compute lower bound transmission t_LB(x)
    I_norm = double(I) ./ reshape(A, 1, 1, 3);
    t_LB = 1 - min(I_norm, [], 3);

    % Impose lower bound on estimated transmission
    t_LB = max(t_LB, 0.1); % impose minimal value of 0.1 and max value of 1
    t_est = min(t_est, 1); % ensure t_est does not exceed 1
    t_est = max(t_est, t_LB);
    
    % Solve optimization problem (Eq. (15))
    % find bin counts for reliability - small bins (#pixels<50) do not comply with 
    % the model assumptions and should be disregarded
    bin_count       = accumarray(idx,1,[clusterCnt,1]);
    bin_count_map   = reshape(bin_count(idx),h,w);
    bin_eval_fun    = @(x) min(1, x/50);

    % Calculate std - this is the data-term weight of Eq. (15)
    K_std = accumarray(idx,R(:),[clusterCnt,1],@std);
    radius_std = reshape( K_std(idx), h, w);
    radius_eval_fun = @(r) min(1, 3*max(0.001, r-0.1));
    radius_reliability = radius_eval_fun(radius_std./max(radius_std(:)));
    data_term_weight   = bin_eval_fun(bin_count_map).*radius_reliability;
    lambda = 0.1;
    t_est = wls_optimization(t_est, data_term_weight, I, lambda);

% Recover dehazed image
t3 = repmat(1.05*t_est, [1 1 3]);
J = I - (1 - t3) .* reshape(A, 1, 1, 3) ./ max(t3, 0.1);

J = min(max(J, 0), 1); % ensure pixel values are in [0, 1]

% contrast enhancement
J = imadjust(J, stretchlim(J), []);

end

function [R, Phi, Theta] = rgb2spherical(IA)
    % Convert RGB to spherical coordinates
    R = sqrt(sum(IA.^2, 3));
    Phi = atan2(IA(:,:,2), IA(:,:,1));
    Theta = acos(IA(:,:,3) ./ (R + eps));
end

function cost = costFun(t_final, t_LB, sigma, lambda, I)

    firstTerm = 0;
    secondTerm = 0;
    % Compute the regularization term (secondTerm)
    [hgt, wdt] = size(t_final);
    for x = 1:hgt
        for y = 1:wdt
            firstTerm = firstTerm + (t_final(x,y) - t_LB(x,y))^2 / (sigma(x,y)^2);
            neighbors = [x-1 y; x+1 y; x y-1; x y+1];
            for n = 1:size(neighbors,1)
                xn = neighbors(n,1);
                yn = neighbors(n,2);
                if xn >= 1 && xn <= hgt && yn >= 1 && yn <= wdt
                    % avoid out-of-bounds
                    diff_t = t_final(x,y) - t_final(xn,yn);
                    diff_I = double(I(x,y,:)) - double(I(xn,yn,:));
                    norm_diff_I = norm(diff_I(:))^2 + eps;
                    secondTerm = secondTerm + (diff_t^2) / norm_diff_I;
                end
            end
        end
    end
    cost = firstTerm + lambda * secondTerm;
end
