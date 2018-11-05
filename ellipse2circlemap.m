% MATLAB verions pre 2018a do not have built-in support Jacobi elliptic 
% functions of complex input arguments -- this code calculates them.
% This code was was taken from: 
% https://au.mathworks.com/matlabcentral/fileexchange/17747-jacobi-elliptic-functions-sn--cn-and-dn-of-complex-phase?focused=5121564&tab=function

%%
function [transform,derivative_transform] = ellipse2circlemap(z,m,a,b,tol)

k = sqrt(m);
if nargin<5, tol = eps; end

% This is the argument of the Jacobi elliptic sine function, it is set up
% to map the ellipse to the circle here, and it will be weighted by the
% sqrt of k once the sn has been calculated.

u =2*ellipke(m)/pi * asin(z./sqrt(a^2 -b^2)); 

% capture memory and save the structure of input arrays
sni = zeros(size(u));
cni = sni;     
dni = sni;

% make a row vector
m = m(:).'; 
u = u(:).';

% represent u in the form u = phi + i*psi
phi = real(u);
psi = imag(u);

[s,c,d] = ellipj(phi,m,tol);
[s1,c1,d1] = ellipj(psi,1-m,tol);

% function evaluations
delta = c1.^2 + m.*s.^2.*s1.^2;
sni(:) = (s.*d1 + sqrt(-1).*c.*d.*s1.*c1)./delta;
cni(:) = (c.*c1 - sqrt(-1).*s.*d.*s1.*d1)./delta;
dni(:) = (d.*c1.*d1 - sqrt(-1).*m.*s.*c.*s1)./delta;

%Check some known identities to make sure things are working as expected:
%These gave zero to within 10^{-15} (as they should)
%max( cni - sqrt(1- sni.^2))
%max( dni - sqrt(1- k.^2*sni.^2))

% weights the function and returns
transform = sqrt(k)*sni;
derivative_transform = sqrt(k)*(2*ellipke(m)/pi./sqrt(a^2-b^2-z.^2)).*dni.*cni;