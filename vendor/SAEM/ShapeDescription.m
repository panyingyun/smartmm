function [SH, coef, ReV] = ShapeDescription(V, S, d)
theta = acos(S(:,3)); 
phi = atan2(S(:,2), S(:,1)); 

SH = RealSphericalHarmonic(theta, phi, d);
coef = (SH'*SH)\(SH'*V);
ReV = SH*coef;
ReV = real(ReV);
end