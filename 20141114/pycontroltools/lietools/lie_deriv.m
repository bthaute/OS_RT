function res = lie_deriv(sf,vf, x)
% lie_deriv(sf,vf, x) calculates the Lie derivative of
% a scalar field sf along a vector field vf
% a vector field :math:`f(x)`

res = jacobian(sf, x)*vf;