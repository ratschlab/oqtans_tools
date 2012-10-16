function score = path_score(state_seq, score_matrix, trans_matrix)

idx = sub2ind(size(score_matrix), state_seq, [1:length(state_seq)]);
score = sum(score_matrix(idx));

idx = sub2ind(size(trans_matrix), state_seq(1:end-1), state_seq(2:end));
score = score + sum(trans_matrix(idx));

