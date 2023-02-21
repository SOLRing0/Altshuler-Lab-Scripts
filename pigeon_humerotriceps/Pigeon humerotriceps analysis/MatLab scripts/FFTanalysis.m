function [FFTdat, max_freq, max_power]=FFTanalysis(Position,Fs,min_PSD)
        L = length(Position); % length of measurement interval
        NFFT = 2^nextpow2(L); % next power of 2 to satisfy the length requirements by either padding with zeros or truncating the signal
        position = Position - mean(Position); % remove DC component and center around zero.
        Y = fft(position,NFFT)/L; % Matlab fft output of the position vector (correction of fft())
        f = Fs/2*linspace(0,1,NFFT/2+1); % creates a vector of reference frequencies and is used for the x-axis
        % %linspace generates evenly spaced points between two specified numbers, in this case 0 and 1 and the third position gives the number of points.
         fourier_function = abs(Y(1:NFFT/2+1)); % amplitude computed 
        % % from determining abs(magnitude), since output is a complex number & 
        % % 2* because output is double-sided, therefore look at one side and
        % % multiply by 2. RMS = voltage amplitude/sqrt(2)
        power = (1/(Fs*L))*abs(fourier_function).^2;
        power(2:end-1) = 2*(power(2:end-1));
        power = power';
        % find peak frequency and power
        max_val = max(2*(fourier_function));
        max_ind = find((2*(fourier_function))>=max_val);
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
FFTdat.power=power;
FFTdat.f=f;
FFTdat.locs=locs;
FFTdat.pks=pks;
FFTdat.ffpow=ffpow;
FFTdat.position=position;


