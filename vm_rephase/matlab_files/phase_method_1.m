function phase_method_1(cptPts, array, params)
%% ESTIMATECOMPLETESIGNAL
% This function estimate the diffuse components and computes the complete
% signal (direct + diffuse).

fprintf('Estimate the full signal (direct + diffuse)...\n');
fLen = params.fLen;
tLen = params.tLen;

arrayCenter = cell2mat(array.center);  % centers
arrayCenter = arrayCenter(:,1:2);
distArrayCptPts = pdist2(cptPts.position,  arrayCenter);
[~, closestArray]= min(distArrayCptPts, [], 2);

mean3 = @(x) sqrt(mean(abs(x).^2,3));
totalDiffuseSTFT = cellfun(mean3, array.meanDiffuseSTFT, 'UniformOutput', false);

omega = 2*pi*params.frequency;
diffuseContribution = zeros(size(cptPts.directEstimate));
for mm = 1:cptPts.N
    fprintf("vm: " + num2str(mm) + "\n");

    % BUILD AMPLITUDE
    distFactor = (1 ./ distArrayCptPts(mm,:).');
    distFactor = distFactor ./ sum(distFactor);
    diffuseSum = zeros(fLen, tLen);

    for dd = 1:array.N
        diffuseSum = diffuseSum + (distFactor(dd) * totalDiffuseSTFT{dd});
    end

    % BUILD PHASE : take the phase of the the closest mic
    arrayIdx = closestArray(mm);
    positions_within_curr_mic_array = cell2mat(array.position(arrayIdx));
    distMics_CptPts = pdist2(cptPts.position(mm,1:2),  positions_within_curr_mic_array(:, 1:2));
    [~, closestMic]= min(distMics_CptPts, [], 2); % must be 1-2-3-4

     % finds the delay between vm and array mic
    pos_closestmic = array.position{arrayIdx}(closestMic, 1:2);
    pos_vm = cptPts.position(mm,1:2);

    d = norm(pos_closestmic - pos_vm);
    mics_vm = exp(1i * d * omega/params.c);  % -?
    mics_vm_f_t = repmat(mics_vm, 1, params.tLen); % flen, tlen

    rephased_stft = angle(array.meanDiffuseSTFT{arrayIdx}(:,:,closestMic) .* mics_vm_f_t);  % flen, tlen

    diffuseSum = diffuseSum .* exp(1i*rephased_stft);

    diffuseContribution(:,mm) = istft(diffuseSum, params.analysisWin, ...
           params.synthesisWin, params.hop, params.Nfft, params.Fs);

end

completeEstimate = cptPts.directEstimate + diffuseContribution;

export_path = 'export/completeEstimate_m1.mat';
save(export_path, 'completeEstimate', '-v7.3'); % time
fprintf("exporting...\n");

end