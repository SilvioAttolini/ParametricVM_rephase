function [completeEstimate, diffuseContribution] = estimateCompleteRephased(cptPts, array, params)
%% ESTIMATECOMPLETESIGNAL
% This function estimate the diffuse components and computes the complete
% signal (direct + diffuse).

fprintf('Estimating the full signals (direct + diffuse)...\n');
fLen = params.fLen;
tLen = params.tLen;

arrayCenter = cell2mat(array.center);  % centers
arrayCenter = arrayCenter(:,1:2);
distArrayCptPts = pdist2(cptPts.position,  arrayCenter);

meanDiffuse = array.meanDiffuse;
for mm = 1:cptPts.N
    fprintf("vm: " + num2str(mm) + "\n");

    pos_vm = cptPts.position(mm,1:2);

    distFactor = (1 ./ distArrayCptPts(mm,:).');
    distFactor = distFactor ./ sum(distFactor);

    diffuseSum = zeros(fLen, tLen);
    for arr = 1:array.N
        for amic = 1:array.micN
            pos_rm = array.position{arr}(amic, 1:2);

            rephasedMeanDiffuse{arr}(:,amic) = getRephased(params, meanDiffuse{arr}(:,amic), pos_rm, pos_vm);
%            display(size(rephasedMeanDiffuse{arr}));

            array.rephasedMeanDiffuseSTFT{arr}(:,:,amic) = stft(rephasedMeanDiffuse{arr}(:,amic), params.analysisWin, ...
                params.hop, params.Nfft, params.Fs);
%            display(size(array.rephasedMeanDiffuseSTFT{arr}(:,:,amic)));

        end
        totalDiffuseSTFT{arr} = sqrt(mean(abs(array.rephasedMeanDiffuseSTFT{arr}).^2,3));
%        display(size(totalDiffuseSTFT{arr}));

        diffuseSum = diffuseSum + (distFactor(arr) * totalDiffuseSTFT{arr});
    end

    diffuseContribution(:,mm) = istft(diffuseSum, params.analysisWin, params.synthesisWin, params.hop, ...
        params.Nfft, params.Fs);

end

completeEstimate = cptPts.directEstimate + diffuseContribution;

export_path = 'export/completeEstimate.mat';
save(export_path, 'completeEstimate', '-v7.3'); % time
fprintf("exporting...\n");

end
