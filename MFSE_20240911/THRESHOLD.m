function [ePhiProj,edproj] = THRESHOLD(ePhi, beta)
%% SIGMOID投影
ePhiProj = 1./(1+exp(-beta*ePhi));
edproj = beta * ePhiProj.*(1-ePhiProj);
end