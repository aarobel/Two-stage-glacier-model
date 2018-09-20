function fout = Base(x,parameters)
%elevation of base below sea level; needs to be consistent with that used
%in evolution problem for h
fout=nan.*ones(size(x));

sill_length = parameters.sill_max-parameters.sill_min;

fout = parameters.icedivide + (parameters.bedslope*x);


fout(x>parameters.sill_min & x<parameters.sill_max) = parameters.icedivide +...
    (parameters.bedslope*parameters.sill_min) +...
    parameters.sill_slope.*(x(x>parameters.sill_min & x<parameters.sill_max)-parameters.sill_min);

fout(x>=parameters.sill_max) = parameters.icedivide +...
    (parameters.bedslope*parameters.sill_min) +...
    parameters.sill_slope.*sill_length + ...
    parameters.bedslope*(x(x>=parameters.sill_max)-parameters.sill_max);

end
