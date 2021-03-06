function fout = dBasedx(x,parameters)
%derivative of elevation of base below sea level; needs to be consistent with that used
%in evolution problem for h

fout = parameters.bedslope.*ones(size(x));
fout(x>parameters.sill_min & x<parameters.sill_max) = parameters.sill_slope;

end
