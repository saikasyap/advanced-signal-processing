
%wvtool(rectwin(48))
%wvtool(bartlett(48))
%wvtool(hamming(48))
%wvtool(hann(48))



w = kaiser(48,4.86);
wvtool(w)


