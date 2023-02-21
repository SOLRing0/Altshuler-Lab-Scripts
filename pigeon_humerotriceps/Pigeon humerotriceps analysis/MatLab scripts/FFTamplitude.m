function [FFTadat]=FFTamplitude(Position,Fs,L,Time)
        NFFT = 2^nextpow2(L); % next power of 2 to satisfy the length requirements by either padding with zeros or truncating the signal
        time = (Time(1:L)-min(Time(1:L)))*1000; % shift start time to 0 and convert to ms
        position = Position(1:L) - mean(Position(1:L)); % remove DC component and center around zero.
        Y = fft(position,NFFT)/L; % Matlab fft output of the position vector (correction of fft())
        f = Fs/2*linspace(0,1,NFFT/2+1); % creates a vector of reference frequencies and is used for the x-axis
        % %linspace generates evenly spaced points between two specified numbers, in this case 0 and 1 and the third position gives the number of points.
         amplitude = abs(Y(1:NFFT/2+1)); % amplitude computed 
        % % from determining abs(magnitude), since output is a complex number & 
        % % 2* because output is double-sided, therefore look at one side and
        % % multiply by 2. RMS = voltage amplitude/sqrt(2)

    %% save parameters
FFTadat.amplitude=amplitude;
FFTadat.f=f;
FFTadat.position=position;
FFTadat.time=time;


