function [GV_cor err alphas] = amp_correct(avgGV, GV, amp_term, alpha_range, alpha_search_grid)
%% Function to find the best amplitude correction term to make the phase velocity from single event similar
%  to average phase velocity
%  written by Ge Jin, jinwar@gmail.com  2013-03-20

alphas = alpha_range(1):alpha_search_grid:alpha_range(2);

for i=1:length(alphas)
	err(i) = ampcorerr(avgGV, GV, amp_term, alphas(i));
end

[temp best_alpha_i] = min(err);

GV_cor=((GV).^-2 + alphas(best_alpha_i).*amp_term').^-.5;

end
