% Region-level overlap analysis of extreme FC deviations (|Z| > thr).
% For each region, count how many edges carry an extreme deviation
% (the region's "deviation degree"). For each degree threshold j = 1..t,
% compute the proportion of subjects with degree >= j, integrate across
% thresholds (AUC), and test the ASD - TD AUC difference against a
% label-shuffling null with FDR.
%
% This demo uses random Z-scores so the script runs standalone. Swap the
% synthetic block for your own normative model output to use on real data.

clear; clc; rng(42);

% Parameters
nRoi  = 390;
K     = nRoi;
thr   = 2.3;
t     = 20;                         % max deviation-degree threshold
nP    = 10000;

% --- Synthetic data (replace with real Z-scores) -------------------------
% In the real pipeline z_asd and z_td come from a GPR normative model
% (PCNtoolkit);
n_asd  = 759;
n_td   = 504;
nEdges = nRoi * (nRoi - 1) / 2;

z_asd = randn(n_asd, nEdges);
z_td  = randn(n_td,  nEdges);
% ------------------------------------------------------------------------

% For each subject and region, count extreme positive / negative
% deviations to that region (nodal degree of the binarised
% deviation adjacency matrix).
[degree_pos_asd, degree_neg_asd] = compute_deviation_degree(z_asd, thr, nRoi);
[degree_pos_td,  degree_neg_td]  = compute_deviation_degree(z_td,  thr, nRoi);

% Stack groups into a single [Nsubs x nRoi] matrix.
% The "- 1" shift lets the j = 1..t threshold loop below evaluate "strictly
% more than one extreme edge" starting from j = 1.
z_idx           = logical([ones(n_asd,1); zeros(n_td,1)]);
degree_pos_orig = [degree_pos_asd - 1, degree_pos_td - 1]';
degree_neg_orig = [degree_neg_asd - 1, degree_neg_td - 1]';
Nsubs           = size(degree_pos_orig, 1);

cnt_neg = zeros(1, K);
cnt_pos = zeros(1, K);

% Permutation test. First iteration is the observed labelling; the rest
% shuffle subjects across the combined stack.
tic
for p = 1:nP
    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end
    asd_mask = z_idx(idx);
    td_mask  = ~asd_mask;

    overlap_pos_asd = zeros(t, K);
    overlap_neg_asd = zeros(t, K);
    overlap_pos_td  = zeros(t, K);
    overlap_neg_td  = zeros(t, K);

    % For each degree threshold j, compute the proportion of subjects per
    % region whose deviation degree meets or exceeds j. Sweeping j gives
    % an overlap-vs-degree curve for each region.
    for j = 1:t
        degree_pos = degree_pos_orig >= j;
        degree_neg = degree_neg_orig >= j;

        overlap_pos_asd(j,:) = sum(degree_pos(asd_mask,:)) / n_asd;
        overlap_neg_asd(j,:) = sum(degree_neg(asd_mask,:)) / n_asd;
        overlap_pos_td(j,:)  = sum(degree_pos(td_mask,:))  / n_td;
        overlap_neg_td(j,:)  = sum(degree_neg(td_mask,:))  / n_td;
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
        % Store the observed statistics and curves for later plotting
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
toc

% Uncorrected permutation p-values, then BH-FDR across regions
punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;
[p_neg]  = mafdr(punc_neg, 'BHFDR', true);
[p_pos]  = mafdr(punc_pos, 'BHFDR', true);


%% Helper function
% For each subject: threshold edges at |Z| > thr, build the binary
% deviation adjacency matrix, and return each region's nodal degree
% (number of extreme-deviation edges attached to it), separately for
% positive and negative deviations.
function [degree_pos, degree_neg] = compute_deviation_degree(z, thr, nRoi)
    z_log_pos = z > thr;
    z_log_neg = z < -thr;

    % Vector -> nRoi x nRoi x nSubs
    adj_pos = recon3dMat(z_log_pos', nRoi);
    adj_neg = recon3dMat(z_log_neg', nRoi);

    nSubs = size(adj_pos, 3);
    degree_pos = zeros(nRoi, nSubs);
    degree_neg = zeros(nRoi, nSubs);

    for i = 1:nSubs
        degree_pos(:,i) = degrees_und(adj_pos(:,:,i));
        degree_neg(:,i) = degrees_und(adj_neg(:,:,i));
    end
end
