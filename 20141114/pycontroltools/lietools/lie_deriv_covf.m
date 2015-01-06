function res = lie_deriv_covf(w,f, x, transpose)
% lie_deriv_covf(w,f, x, transpose) calculates the Lie derivative of the
% covector field w along the vector field f 
% includes the option to omit the transposition of (jacobian(w.' ,x))

if transpose == true
    res = f.' * (jacobian(w.' ,x)).' + w * jacobian(f,x);
end
    
if transpose == false
    res = f.' * (jacobian(w.' ,x)) + w * jacobian(f,x);
end