function [sig_out] = Shuffle_ISI(sig_in)

if size(sig_in,1) > size(sig_in,2)
    sig_in = sig_in';
end
svec = sig_in>0;
if any(svec)
    rstart = randi(find(svec,1,'first')-1);
    dsig_start = [rstart find(svec)];
    dsig_end = [find(svec)-1];
    if length(dsig_start) ~= length(dsig_end)
        dsig_start  = dsig_start(1:min([length(dsig_start) length(dsig_end)]));
        dsig_end    = dsig_end(1:min([length(dsig_start) length(dsig_end)]));
    end
    sig_out = zeros(length(svec),1); % randomly shuffle interspike intervals
    [~, rand_order] = sort(rand(length(dsig_start),1));
    dsig_start = dsig_start(rand_order);
    dsig_end = dsig_end(rand_order);
    rand_inds = cumsum([dsig_end+1 - dsig_start]);
    sig_out(rand_inds) = 1;
else
    sig_out = sig_in;
end

end