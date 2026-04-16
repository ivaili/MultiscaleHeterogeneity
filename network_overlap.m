% Network-level overlap analysis of extreme FC deviations (|Z| > thr).
% Edges are grouped into within- and between-network blocks defined by a
% 10-network partition (7 cortical Yeo + thalamus, striatum, amyg/hipp),
% giving 55 unique blocks. For each block and each "network-deviation"
% threshold, we compute the proportion of subjects whose block deviation
% exceeds that threshold, integrate over thresholds (AUC), and test the
% ASD - TD AUC difference against a label-shuffling null with FDR.
%
% This demo uses random Z-scores so the script runs standalone. Swap the
% synthetic block for your own normative model output to use on real data.

clear; clc; rng(42);

% Parameters
nRoi       = 390;
K          = 55;                    % 10 networks -> 55 unique blocks
thr        = 2.3;
nP         = 10000;
thresholds = 0.01:0.01:0.1;         % block-level deviation thresholds
nThr       = length(thresholds);

% --- Synthetic data (replace with real Z-scores) -------------------------
% In the real pipeline z_asd and z_td come from a GPR normative model
% (PCNtoolkit)
n_asd = 759;
n_td  = 504;
nEdges = nRoi * (nRoi - 1) / 2;

z_asd = randn(n_asd, nEdges);
z_td  = randn(n_td,  nEdges);

% Generate 10-network partition: split the 390 regions roughly evenly.
% Replace with the actual `networks` vector for analyses.
networks = repelem(1:10, ceil(nRoi/10));
networks = networks(1:nRoi)';
netNames = {'Visual','Somatomotor','Dorsalattn','Ventralattn','Limbic', ...
            'Frontoparietal','Defaultmode','Thalamus','Striatum','AmygHypp'};
% ------------------------------------------------------------------------

% For each subject, aggregate edge-level extreme deviations into the 55
% within/between-network blocks (separately for positive and negative).
[z_pos_asd, z_neg_asd] = classify_edges(z_asd, thr, networks, netNames);
[z_pos_td,  z_neg_td]  = classify_edges(z_td,  thr, networks, netNames);

% Stack groups for the permutation
z_idx      = logical([ones(n_asd,1); zeros(n_td,1)]);
z_all_pos  = [z_pos_asd; z_pos_td];
z_all_neg  = [z_neg_asd; z_neg_td];
Nsubs      = size(z_all_pos,1);

cnt_neg = zeros(1, K);
cnt_pos = zeros(1, K);

% Permutation test. First iteration is the observed labelling; the rest
% shuffle subjects across the combined stack.
for p = 1:nP

    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end

    asd_mask = z_idx(idx);
    td_mask  = ~asd_mask;

    overlap_pos_asd = zeros(nThr, K);
    overlap_neg_asd = zeros(nThr, K);
    overlap_pos_td  = zeros(nThr, K);
    overlap_neg_td  = zeros(nThr, K);

    % For each block-deviation threshold, compute the proportion of
    % subjects whose block-level deviation score meets or exceeds it.
    % Sweeping across thresholds gives an overlap-vs-threshold curve.
    for counter = 1:nThr
        j = thresholds(counter);
        netdev_pos = z_all_pos >= j;
        netdev_neg = z_all_neg >= j;

        overlap_pos_asd(counter,:) = sum(netdev_pos(asd_mask,:)) / n_asd;
        overlap_neg_asd(counter,:) = sum(netdev_neg(asd_mask,:)) / n_asd;
        overlap_pos_td(counter,:)  = sum(netdev_pos(td_mask,:))  / n_td;
        overlap_neg_td(counter,:)  = sum(netdev_neg(td_mask,:))  / n_td;
    end

    % Summarise each curve by its AUC (trapezoidal), then take the ASD-TD
    % difference as the test statistic.
    auc_pos_asd = trapz(overlap_pos_asd);
    auc_neg_asd = trapz(overlap_neg_asd);
    auc_pos_td  = trapz(overlap_pos_td);
    auc_neg_td  = trapz(overlap_neg_td);

    delta_auc_pos = auc_pos_asd - auc_pos_td;
    delta_auc_neg = auc_neg_asd - auc_neg_td;

    if p == 1
        % Store the observed statistics and the curves for later plotting
        delta_auc_pos1   = delta_auc_pos;
        delta_auc_neg1   = delta_auc_neg;
        perc_diff_pos    = (delta_auc_pos ./ auc_pos_td) * 100;
        perc_diff_neg    = (delta_auc_neg ./ auc_neg_td) * 100;
        overlap_pos_asd1 = overlap_pos_asd;
        overlap_neg_asd1 = overlap_neg_asd;
        overlap_pos_td1  = overlap_pos_td;
        overlap_neg_td1  = overlap_neg_td;
    else
        cnt_neg = cnt_neg + (delta_auc_neg >= delta_auc_neg1);
        cnt_pos = cnt_pos + (delta_auc_pos >= delta_auc_pos1);
    end

end

% Uncorrected permutation p-values, then BH-FDR across the 55 blocks
punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;
[p_neg]  = mafdr(punc_neg, 'BHFDR', true);
[p_pos]  = mafdr(punc_pos, 'BHFDR', true);


%% Helper function
% Takes a subjects-by-edges matrix of Z-scores, thresholds into extreme
% positive/negative deviations, maps them onto the 10 x 10 network grid
% via plotClassifiedEdges2, and returns the flattened 55-block vector
% per subject.
function [z_pos, z_neg] = classify_edges(z, thr, networks, netNames)
    nRoi = 390;
    z_log_pos = z > thr;
    z_log_neg = z < -thr;

    % Vector -> 390 x 390 x nSubs
    z_log_pos_sq = recon3dMat(z_log_pos', nRoi);
    z_log_neg_sq = recon3dMat(z_log_neg', nRoi);

    nSubs = size(z_log_pos_sq, 3);
    z_net_pos = zeros(10, 10, nSubs);
    z_net_neg = zeros(10, 10, nSubs);

    for i = 1:nSubs
        [~,~,z_net_pos(:,:,i)] = plotClassifiedEdges2(z_log_pos_sq(:,:,i), networks, 0, netNames);
        [~,~,z_net_neg(:,:,i)] = plotClassifiedEdges2(z_log_neg_sq(:,:,i), networks, 0, netNames);
    end

    % 10 x 10 x nSubs -> nSubs x 55 (lower triangle incl. diagonal)
    z_pos = flatten3dMat(z_net_pos, 0)';
    z_neg = flatten3dMat(z_net_neg, 0)';
end
