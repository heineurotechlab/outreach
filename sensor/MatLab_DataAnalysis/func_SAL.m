
function [SAL] = func_SAL(vx,vy,vz)
    speed = sqrt(vx.^2 + vy.^2+vz.^2);
    Ts = 1/1000; % sample time, seconds
    %% Bits taken from script of Balasubramanian et al. 2012
    parameters = [20,4]; % Recommended in SAL paper, Balasubramanian et al. 2012
    % Calculate the spectrum of the speed profile.
    N = length(speed);
    Nfft = 2^(ceil(log2(N))+parameters(2));
    speedSpectrum = abs(fft( speed, Nfft ));

    % Normalize spectrum with respect to the DC component.
    speedSpectrum = speedSpectrum/speedSpectrum(1);
    % Get index corresponding to the cut off frequency.
    freq = 0:(1/Ts)*(1/Nfft):(1/Ts)*((Nfft-1)/Nfft);
    inxFc = find( freq(1:end) <= parameters(1), 1, 'last' );

    % Calculate the spectral arc length.
    % 1. select the spectrum of interest.
    speedSpectrum = speedSpectrum(1:inxFc);
    % 2. Calculate the incremental arc lengths.
    dArcLengths = sqrt((1/(inxFc-1))^2 + (diff(speedSpectrum)).^2);
    % 4. Compute movement smoothness.
    SAL = -sum(dArcLengths);
end