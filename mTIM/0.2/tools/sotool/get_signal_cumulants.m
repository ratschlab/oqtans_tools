function [m,s,b] = get_signal_cumulants(sig)

DIM = size(sig,1);
m = zeros(DIM,1);
s = zeros(DIM,1);
b = zeros(DIM,1);

for d=1:DIM,
    idx = find(~isinf(sig(d,:)) & ~isnan(sig(d,:)));
    if ~isempty(idx),
        m(d) = mean(sig(d,idx));
        s(d) = std(sig(d,idx));
        b(d) = 1; % valid signal
    end
end
