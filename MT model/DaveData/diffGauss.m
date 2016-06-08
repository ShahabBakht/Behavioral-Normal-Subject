function ri=diffGauss(a,data)
% a(1) = excitation amplitude
% a(2) = excitation size
% a(3) = inhibition amplitude
% a(4) = inhibition size
% a(5) = baseline response

%ri = a(1)*erf(data/a(2))-a(3)*erf(data/(a(4)))+a(5);
ri = a(1)*erf(data/a(2))-a(3)*erf(data/(a(2)+a(4)))+a(5);