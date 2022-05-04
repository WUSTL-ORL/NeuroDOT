function Rmap=FisherZ2R(Zmap)

% This function transforms z-value data into r-value data via the
% transform: r = [exp(2z)-1]/[exp(2z)+1];


Rmap=((exp(2.*Zmap))-1)./((exp(2.*Zmap))+1);
Rmap(~isfinite(Rmap))=0;