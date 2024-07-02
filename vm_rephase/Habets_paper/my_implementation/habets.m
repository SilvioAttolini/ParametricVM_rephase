close all;
clear;

%   import array, cptPts, params
load('data/cptPts.mat', 'cptPts');
load('data/array.mat', 'array');
load('data/params.mat', 'params');

display(params);
display(cptPts);
display(array);

% Initialization
%Fs = 8000; % Sample frequency (Hz)
%c = 340; % Sound velocity (m/s)
%K = 256; % FFT length
%M = 4; % Number of sensors
%d = 0.1; % Inter sensor distance (m)
%type_nf = 'spherical'; % Type of noise field: 'spherical' or 'cylindrical'
%L = 20*Fs; % Data length

Fs = params.Fs; % /2;  % for now, then remove /2
c = params.c;
K = 256;
M = 2; % calculate the SC between 1 vm couple at a time cptPts.N;
%d = norm(cptPts.position(1, :) - cptPts.position(2, :));
d = 0.2;
type_nf = 'spherical'; % Type of noise field
L = size(cptPts.directEstimate);
L = L(1);

display(d);
display(L);

%% Generate M mutually 'independent' babble speech input signals
%[data,Fs_data] = audioread('babble_8kHz.wav');
%if Fs ~= Fs_data
%    error('Sample frequency of input file is incorrect.');
%end
% instead of reading from a was, use directly the noise = complete - direct of the arrays

%mean3 = @(x) sqrt(mean(abs(x).^2,3));
%totalDiffuseSTFT = cellfun(mean3, array.meanDiffuseSTFT, 'UniformOutput', false);
%totalDiffuse = istft(totalDiffuseSTFT{1});
%data = totalDiffuse;
data1 = array.meanDiffuse{1};
data1 = data1(:, 1);
data5 = array.meanDiffuse{5};
data5 = data5(:, 1);

data1 = data1(:); % Ensure data1 is a column vector
data5 = data5(:); % Ensure data2 is a column vector
data = [data1; data5];
display(size(data));

% Combine data1 and data5 into a stereo signal
stereoData = [data1, data5];
fileName = 'stereo.wav';
audiowrite(fileName, stereoData, Fs);
disp(['WAV file created: ', fileName]);

data = data - mean(data);

% data = mean(data, 2);
%data = [data; data]; nope! creates SC = 1 because of the diagonal=1?


babble = zeros(L,M);
for m=1:M
    babble(:,m) = data((m-1)*L+1:m*L);
end

%% Generate matrix with desired spatial coherence
ww = 2*pi*Fs*(0:K/2)/K;
DC = zeros(M,M,K/2+1);
for p = 1:M
    for q = 1:M
        if p == q
            DC(p,q,:) = ones(1,1,K/2+1);
        else
            % case 'spherical'
            DC(p,q,:) = sinc(ww*abs(p-q)*d/(c*pi));
        end
    end
end

%% Generate sensor signals with desired spatial coherence
%x = mix_signals(babble,DC,'cholesky');
x = mix_signals(babble,DC,'eigen');

%% Compare desired and generated coherence
K_eval = 256;
ww = 2*pi*Fs*(0:K_eval/2)/K_eval;
sc_theory = zeros(M-1,K_eval/2+1);
sc_generated = zeros(M-1,K_eval/2+1);

% Calculalte STFT and PSD of all output signals
X = stft(x,'Window',hanning(K_eval),'OverlapLength',0.75*K_eval,'FFTLength',K_eval,'Centered',false);
X = X(1:K_eval/2+1,:,:);
phi_x = mean(abs(X).^2,2);

% Calculate spatial coherence of desired and generated signals
for m = 1:M-1
    % case 'spherical'
    sc_theory(m,:) = sinc(ww*m*d/(c*pi));

    % Compute cross-PSD of x_1 and x_(m+1)
    % psi_x =  mean(X(:,:,1) .* conj(X(:,:,m+1)),2);
    % Compute cross-PSD of x_m and x_(m+1), odd-even
    psi_x =  mean(X(:,:,1) .* conj(X(:,:,m+1)),2);

    % Compute real-part of complex coherence between x_1 and x_(m+1)
    sc_generated(m,:) = real(psi_x ./ sqrt(phi_x(:,1,1) .* phi_x(:,1,m+1))).';
end

% Calculate normalized mean square error
NMSE = zeros(M,1);
for m = 1:M-1
    NMSE(m) = 10*log10(sum(((sc_theory(m,:))-(sc_generated(m,:))).^2)./sum((sc_theory(m,:)).^2));
end

% Plot spatial coherence of the pair
figure(1);
% MM=min(2,M-1);
MM = M-1;
Freqs=0:(Fs/2)/(K/2):Fs/2;
for m = 1:MM
    subplot(MM,1,m);
    plot(Freqs/1000,sc_theory(m,:),'-k','LineWidth',1.5)
    hold on;
    plot(Freqs/1000,sc_generated(m,:),'-.b','LineWidth',1.5)
    hold off;
    xlabel('Frequency [kHz]');
    ylabel('Real(Spatial Coherence)');
    title(sprintf('Inter sensor distance %1.2f m',m*d));
    legend('Theory',sprintf('Proposed Method (NMSE = %2.1f dB)',NMSE(m)));
    grid on;
    saveas(gcf, 'habets_SC_eigen.png');
end

% Save babble speech
audiowrite('habets_result.wav',x,Fs);
%
%pause(99);
