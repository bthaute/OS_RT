function res = lie_bracket(f,g, x)
% lie_bracket(f,g, x) calculates the Lie bracket for the
% vector field g along the vector field f
% along the vector field :math:`f(x)`


res = jacobian(g ,x)*f - jacobian(f,x)*g;