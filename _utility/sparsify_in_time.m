function [sig_out] = sparsify_in_time(t_vec, sig_in, time_win);

sig_in = shockVec;
t_vec = ms.timestamps./1000;
time_win = 5;

