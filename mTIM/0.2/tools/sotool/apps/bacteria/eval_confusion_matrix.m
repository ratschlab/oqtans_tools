
function ret = eval_confusion_matrix(tp, tn, fp, fn)

    sensitivity = tp / (tp + fn); %same as recall
    specificity = tn / (tn + fp);
    precision = tp / (tp + fp); %also called positive predictive value

    % harmonic mean of precision and recall
    f_score = 2 * precision * sensitivity / (precision + sensitivity);

    ret = {};
    ret.sensitivity = sensitivity;
    ret.specificity = specificity;
    ret.precision = precision;
    ret.f_score = f_score;

