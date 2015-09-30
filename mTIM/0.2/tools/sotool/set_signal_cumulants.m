function signal = set_signal_cumulants(signal, m,s,b)


DIM = size(signal,1);

for d=1:DIM,
    idx = find(~isinf(signal(d,:)) & ~isnan(signal(d,:)));
    if ~isempty(idx) && b(d),
        % if signal contains a valid sequence and
        % a conversion of this dimension is defined

        % center signal
        signal(d,idx) = signal(d,idx)-mean(signal(d,idx));
        % default variance
        % check if there is any variance
        ds = std(signal(d,idx));
        if abs(ds)>1e-8,
            signal(d,idx) = signal(d,idx)/ds;
        end       
        % target mean and variance
        signal(d,idx) = signal(d,idx)*s(d);
        signal(d,idx) = signal(d,idx)+m(d);
    end
end


