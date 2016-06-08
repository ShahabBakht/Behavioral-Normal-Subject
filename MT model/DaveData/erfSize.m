function ri=erfSize(a,data)
% a(1) = excitation amplitude
% a(2) = excitation size
% a(3) = baseline response

ri = a(1)*erf(data/a(2))+a(3);