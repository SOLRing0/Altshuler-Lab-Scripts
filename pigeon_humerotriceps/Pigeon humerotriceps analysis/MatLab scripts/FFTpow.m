function [FFTpdat, max_freq, max_power]=FFTpow(Fs,min_PSD,L,f,mean_amp)
        power = (1/(Fs*L))*abs(mean_amp).^2; % convert mean amplitude to PSD
        power(2:end-1) = 2*(power(2:end-1));
        power = power';
        % find peak frequency and power
        max_val = max(2*(mean_amp));
        max_ind = find((2*(mean_amp))>=max_val);
        max_freq = f(max_ind);
        max_power = power(max_ind);
% find all frequencies with power above the minimum threshold
% get locations of PSD values above min and remove any that correspond to
% frequencies near zero (DC)
locs = find(power>=min_PSD);
% find fundamental/contributing frequencies
pks = f(locs);
% find their associated PSD values
ffpow = power(locs);

    %% save parameters
FFTpdat.power=power;
FFTpdat.f=f;
FFTpdat.locs=locs;
FFTpdat.pks=pks;
FFTpdat.ffpow=ffpow;



