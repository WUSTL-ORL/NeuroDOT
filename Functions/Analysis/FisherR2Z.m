function Zmap=FisherR2Z(Rmap)

% This function transforms r-value data into z-value data via the
% transform: z = 1/2 ln([1+r]/[1-r]) = arctanh(r)

%Zmap=arctanh(Rmap);

Zmap=0.5.*(log(1+Rmap)-log(1-Rmap));
Zmap(~isfinite(Zmap))=0;