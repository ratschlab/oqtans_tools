
function [tp, tn, fp, fn] = eval_performance(true_seq, pred_seq, target_state)
    % target state can either be exon, start or stop

    % create a mask based on target state
    true_exon = true_seq == target_state;
    pred_exon = pred_seq == target_state;

    tp = 0;
    tn = 0;
    fp = 0;
    fn = 0;

    for j=1:length(true_seq),

        true_n = true_exon(j);
        pred_n = pred_exon(j);
    
        % fill confusion matrix
        if pred_n == 1 && true_n == 1
            tp = tp + 1;
        elseif pred_n == 0 && true_n == 0
            tn = tn + 1;
        elseif pred_n == 1 && true_n == 0
            fp = fp + 1;
        elseif pred_n == 0 && true_n == 1
            fn = fn + 1;
        end

    end

    
