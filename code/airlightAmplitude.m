function A = estimateAirlightAmplitude(img, Ahat)
% estimateAirlightAmplitude  Recover airlight magnitude a = ‖A‖
% Inputs:
%   img  - H×W×3 hazy RGB image (double in [0,1])
%   Ahat - unit-vector airlight orientation A/‖A‖ (3×1)
% Output:
%   A    - full atmospheric-light vector = a_est * Ahat where Ahat is the
%   airlight orientation



  %--- Step 1: Dehaze assuming a = 1 (init J, t) ------------------------
  % (prepares l_a and t_a before amplitude recovery)
  [J, t] = dehazeHazeLines(img, Ahat); % J is the haze-free image, estimate of the true scene color at x after "removing" airlight



  %--- Step 2: Compute per-pixel brightness l_a(x) (Eq. 6) ------------
  %   l_a(x) = l(x) = ‖J(x)‖₂, with σ_a in Eq. (6) 
  l = sqrt(sum(J.^2, 3));

  %--- Step 3: Build l*_a(s) = percentile_{99%}{l_a | t_a=s} (Eq. 8) --
  numBins = 100;
  sVals   = linspace(0, 1, numBins+1);
  lstar_a = nan(1, numBins);
  halfBin = 1/numBins/2;
  for i = 1:numBins
    mask = abs(t - sVals(i+1)) < halfBin;
    if any(mask(:))
      lstar_a(i) = prctile(l(mask), 99);
    end
  end

  %--- Step 4: Define cost(a,k) implementing Eqs. 9–10 --------------
  costFun = @(ak) amplitudeCost(ak, sVals(2:end), lstar_a);
  %   where Eq. 9 maps s→(s-1)/a+1, and Eq. 10 is the sum of squares

  %--- Step 5: Solve for [a_est, k_est] via fminsearch (Eq. 10), minimize the cost -------
  initParams = [1, 1];
  params = fminsearch(costFun, initParams);
  a_est  = params(1);

  %--- Step 6: Ensure no invalid mapped s ((s-1)/a+1<0) -------------
  s = sVals(2:end);
  while any((s - 1)/a_est + 1 < 0)
    Ahat = 1.2 * Ahat;                 % bump orientation norm
    [J, t] = dehazeHazeLines(img, Ahat);
    l = sqrt(sum(J.^2, 3));
    for i = 1:numBins
      mask = abs(t - sVals(i+1)) < halfBin;
      if any(mask(:))
        lstar_a(i) = prctile(l(mask), 99);
      end
    end
    costFun = @(ak) amplitudeCost(ak, s, lstar_a);
    params  = fminsearch(costFun, params);
    a_est   = params(1);
  end

  %--- Step 7: Form final airlight vector ----------------------------
  A = a_est * Ahat;
end


function cost = amplitudeCost(ak, s, lstar_a)
% amplitudeCost  Compute objective from Eqs. 9–10
%   ak      - [a; k]
%   s       - transmission bins
%   lstar_a - l*_a(s) values

  a = ak(1);  % candidate amplitude
  k = ak(2);  % approx. constant true brightness

  % Eq. 9: map original s → s_a = (s-1)/a + 1
  s_a = (s - 1)/a + 1;

  % interpolate l*_a at s_a
  l_interp = interp1(s, lstar_a, s_a, 'pchip', NaN);

  % Eq. 6: σ_a(s) = a·s/(s + a – 1)
  sigma = a .* s ./ (s + a - 1);

  % Eq. 10: sum of squared errors
  err  = (l_interp - sigma * k).^2;
  cost = nansum(err);
end
