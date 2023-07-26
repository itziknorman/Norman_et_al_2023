function [f_out_log,f_log_ind] = convert_linfreq_to_logfreq(f_out,nfreqs)
tmpfreqs = linspace(log(f_out(1)), log(f_out(end)), nfreqs);
f_out_log_tmp = exp(tmpfreqs);
f_out_log = []; f_log_ind = [];
for k = 1:nfreqs
    [~, f_log_ind(k)] = min(abs(f_out-f_out_log_tmp(k)));
    f_out_log(k) = f_out(f_log_ind(k));
end