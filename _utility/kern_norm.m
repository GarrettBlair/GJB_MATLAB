function dat = kern_norm(dat_in, kern)
kp = [zeros(length(dat_in),1); kern; zeros(length(dat_in),1)];
start = (length(kp)+1)/2;
for i = 1:length(dat_in);
    tempkern = kp(start:start+length(dat_in)-1);
    tempkern = tempkern/sum(tempkern);

    dat(i) = nansum(tempkern'.*dat_in);
    start = start-1;
end
