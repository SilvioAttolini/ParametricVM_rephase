function quickload()

load('data/cptPts.mat', 'cptPts');
load('data/array.mat', 'array');
load('data/params.mat', 'params');

[completeEstimate, diffuseContribution] = estimateCompleteRephased(cptPts, array, params);

end
