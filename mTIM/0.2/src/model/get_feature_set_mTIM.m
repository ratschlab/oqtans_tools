function [FEAT, MAP] = get_feature_set_mTIM()

% FEAT = get_feature_set_mTIM()
%
% Returns a cell with fields specifying the names and order 
% of the features used or learning.
%
% written by Georg Zeller, MPI Tuebingen, Germany, 2008-2011

% W and C denote the forward (Watson) and reverse (Crick) strand 
% of the genome, respectively

FEAT = {'exon_cover', ...       %1
        'exon_diff', ...        %2
        'cover_grad', ...       %3
        'intron_span', ...      %4
        'intron_diff', ...      %5
        'intron_start_W', ...   %6
        'intron_end_W', ...     %7
        'intron_start_C', ...   %8
        'intron_end_C', ...     %9
        'don_pred_W', ...       %10
        'acc_pred_W', ...       %11
        'don_pred_C', ...       %12
        'acc_pred_C', ...       %13
        'pair_span', ...        %14
        'total_cover', ...      %15
        'low_cover_blocks', ... %16
        'repeats' ...           %17
        'intron_span_low' ...   %18
        'intron_span_med' ...   %19
        'intron_span_high' ...  %20
        'cufflinks_feat_W' ...  %21
        'cufflinks_feat_C' ...  %22
        };

% define feature groups 
MAP = {'use_filtered_intron_feats',[4]; ...
       'use_pair_feats',[14:16]; ...
       'use_splice_feats',[10:13]; ...   
       'use_repeat_feats',[17]; ...
       'use_binned_span_feats',[18:20]; ...
       'use_cuffl_feats',[21,22] ...
      };


% eof
