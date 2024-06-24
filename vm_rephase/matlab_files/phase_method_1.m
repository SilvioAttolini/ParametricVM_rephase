function phase_method_1(cptPts, array, params)
%% ESTIMATECOMPLETESIGNAL
% This function estimate the diffuse components and computes the complete
% signal (direct + diffuse).

fprintf('Estimate the full signal (direct + diffuse)...\n');
%fLen = params.fLen;
tLen = params.tLen;

arrayCenter = cell2mat(array.center);  % centers
arrayCenter = arrayCenter(:,1:2);
distArrayCptPts = pdist2(cptPts.position,  arrayCenter);
[~, closestArray]= min(distArrayCptPts, [], 2);

mean3 = @(x) sqrt(mean(abs(x).^2,3));
totalDiffuseSTFT = cellfun(mean3, array.meanDiffuseSTFT, 'UniformOutput', false);

d = 0.078520;

% USING STFT
% omega = 2*pi*params.frequency;

diffuseContribution = zeros(size(cptPts.directEstimate));
for mm = 1:cptPts.N
    fprintf("vm: " + num2str(mm) + "\n");

    % BUILD AMPLITUDE
    distFactor = (1 ./ distArrayCptPts(mm,:).');
    distFactor = distFactor ./ sum(distFactor);
    diffuseSum = zeros(params.fLen, tLen);

%    if mod(mm, 2) == 0
%        % totDiffRephased = totDiff_FFT .* exp(-1i * (distArrayCptPts(mm,dd) + d) .* omega / params.c);
%        % serve una porzione di d in base all'orientamento
%        k = 2 * pi * d / params.c;
%        sinc = sin(k * omega_old) ./ (k * omega_old);
%        display(size(sinc));
%        totDiffRephased_even = totDiffRephased_odd .* exp(-1i .* d);
%
%        diffTotRephased_time = real(ifft(totDiffRephased_even, fLen));
%        diffTotReph_STFT = stft(diffTotRephased_time, params.synthesisWin, params.hop, params.Nfft, params.Fs);
%
%        for dd = 1:array.N
%            diffuseSum = diffuseSum + (distFactor(dd) * diffTotReph_STFT);
%        end
%    else
    for dd = 1:array.N

        totDiff_time = real(istft(totalDiffuseSTFT{dd}, params.analysisWin, ...
            params.synthesisWin, params.hop, params.Nfft, params.Fs));

        % USING FFT
        % fLen = params.Fs * length(diffuseSum_time);
        fLen = length(totDiff_time);
        freqs = linspace(0, params.Fs, fLen);
        % freqs_flips = [freqs(1: end/2), flip(freqs(1: end/2))];
        omega = 2*pi*freqs;
        omega_old = omega;

        totDiff_FFT = fft(totDiff_time, fLen);

        totDiffRephased = totDiff_FFT .* exp(-1i * distArrayCptPts(mm,dd) .* omega / params.c);
        totDiffRephased_odd = totDiffRephased;

        diffTotRephased_time = real(ifft(totDiffRephased, fLen));
        diffTotReph_STFT = stft(diffTotRephased_time, params.synthesisWin, params.hop, params.Nfft, params.Fs);

        diffuseSum = diffuseSum + (distFactor(dd) * diffTotReph_STFT);
       % diffuseSum = diffuseSum + (distFactor(dd) * totalDiffuseSTFT{dd});
    end

    diffuseSum_time = real(istft(diffuseSum, params.analysisWin, ...
        params.synthesisWin, params.hop, params.Nfft, params.Fs));


%    diffuseSum_FFT = fft(diffuseSum_time, fLen);


    % BUILD PHASE : take the phase of the the closest mic
%    arrayIdx = closestArray(mm);
%    positions_within_curr_mic_array = cell2mat(array.position(arrayIdx));
%    distMics_CptPts = pdist2(cptPts.position(mm,1:2),  positions_within_curr_mic_array(:, 1:2));
%    [~, closestMic]= min(distMics_CptPts, [], 2); % must be 1-2-3-4
%
%     % finds the delay between vm and array mic
%    pos_closestmic = array.position{arrayIdx}(closestMic, 1:2);
%    pos_vm = cptPts.position(mm,1:2);
%    d = norm(pos_closestmic - pos_vm);
%    vm_delay = exp(1i * d * omega/params.c);  % -?  % 16000 freqs
    % angle_delay = angle(vm_delay);

    % USING STFT  limited by time-space resolution
%    mics_vm_f_t = repmat(mics_vm, 1, params.tLen); % flen, tlen
%    rephased_stft = angle(array.meanDiffuseSTFT{arrayIdx}(:,:,closestMic) .* mics_vm_f_t);  % flen, tlen
%    diffuseSum = diffuseSum .* exp(1i*rephased_stft);
%    diffuseContribution(:,mm) = istft(diffuseSum, params.analysisWin, ...
%           params.synthesisWin, params.hop, params.Nfft, params.Fs);

    % USING FFT
%    meanDiff_time = real(istft(array.meanDiffuseSTFT{arrayIdx}(:,:,closestMic), params.analysisWin, ...
%           params.synthesisWin, params.hop, params.Nfft, params.Fs));
%    meanDiff_FFT = fft(meanDiff_time, fLen);  % 16000
%    % angles_fft = angle(meanDiff_FFT);
%
%     rephased_fft = angles_fft + angle_delay;  % exp arguments
%    diffuseSum_FFT_phased = diffuseSum_FFT .* vm_delay;
%     diffuseSum_FFT_phased = diffuseSum_FFT .* exp(1i*rephased_fft);

%    diffuseContribution(:, mm) = real(ifft(diffuseSum_FFT_phased, fLen)); % time
    diffuseContribution(:, mm) = diffuseSum_time;
end

completeEstimate = cptPts.directEstimate + diffuseContribution;

export_path = 'export/completeEstimate_prephased.mat';
save(export_path, 'completeEstimate', '-v7.3'); % time
fprintf("exporting...\n");

end