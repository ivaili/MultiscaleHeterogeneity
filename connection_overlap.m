% Edge-level overlap analysis of extreme FC deviations (|Z| > thr).
% For each edge, compare the proportion of ASD vs TD subjects with an
% extreme deviation, and test against a label-shuffling null with FDR.
%
% This demo uses random Z-scores so the script runs standalone. Swap the
% synthetic block for your own normative model output to use on real data.

clear; clc; rng(42);

% Parameters
nRoi  = 390;
K     = nRoi * (nRoi - 1) / 2;   % 75,855 unique edges
thr   = 2.3;
nP    = 10000;
n_asd = 759;
n_td  = 504;

% --- Synthetic data (replace with real Z-scores) -------------------------
% In the real pipeline z_asd and z_td come from a GPR normative model
% (PCNtoolkit) and z_idx is loaded from z_all.mat.
z_asd = randn(n_asd, K);
z_td  = randn(n_td,  K);
z_idx = [true(n_asd, 1); false(n_td, 1)];
% ------------------------------------------------------------------------

z_idx = logical(z_idx);
z_all = [z_asd; z_td];

% Binarise positive and negative extreme deviations separately
z_all_pos = z_all >  thr;
z_all_neg = z_all < -thr;

Nsubs   = size(z_all_pos, 1);
cnt_neg = zeros(1, K);
cnt_pos = zeros(1, K);

% Permutation test. First iteration is the observed (identity) labelling;
% the rest shuffle subjects across the combined stack.
tic
for p = 1:nP
    if p == 1
        idx = 1:Nsubs;
    else
        idx = randperm(Nsubs);
    end
    asd_mask = z_idx(idx);
    td_mask  = ~asd_mask;

    % TD proportion minus ASD proportion, per edge. Scaling by group size
    % matters because the groups are unbalanced.
    overlap_neg = sum(z_all_neg(td_mask,:))  / n_td ...
                - sum(z_all_neg(asd_mask,:)) / n_asd;
    overlap_pos = sum(z_all_pos(td_mask,:))  / n_td ...
                - sum(z_all_pos(asd_mask,:)) / n_asd;

    if p == 1
        overlap1_neg = overlap_neg;
        overlap1_pos = overlap_pos;
    else
        cnt_neg = cnt_neg + (overlap_neg >= overlap1_neg);
        cnt_pos = cnt_pos + (overlap_pos >= overlap1_pos);
    end
end
toc

% Uncorrected permutation p-values, then BH-FDR across all edges
punc_neg = cnt_neg / nP;
punc_pos = cnt_pos / nP;

p_neg = mafdr(punc_neg, 'BHFDR', true);
p_pos = mafdr(punc_pos, 'BHFDR', true);

% Back to matrix form for plotting / network aggregation
p_neg_recon = reconMat(p_neg', nRoi);
p_pos_recon = reconMat(p_pos', nRoi);
