% Function to calculate the difference between the averaged phase velocity and single event phase velocity map
% written by Ge Jin
% jinwar@gmail.com
% Mar, 2012
function err=ampcorerr(avgGV, GV, amp_term, alpha)

	GV_cor=((GV).^-2 + alpha.*amp_term').^-.5;
	
	err = nanmean(nanmean(abs(GV_cor-avgGV)));
		
end


