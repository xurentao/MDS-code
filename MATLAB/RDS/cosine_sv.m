function [cs,sn] = cosine_sv(sv)
% Computes the cosine and sine of the singular values.
%
% Software based on: "Swapping 2x2 blocks in the Schur and generalized
%                     Schur form" - D. Camps, N. Mastronardi, R. Vandebril,
%                     and P. Van Dooren.
k = length(sv);
cs = zeros(k,1);
sn = zeros(k,1);
for i = 1:k
  if sv(i) > 1
    cs(i) = 1/sqrt(1+sv(i)^2);
    sn(i) = cs(i)*sv(i);
  else
    sn(i) = 1/sqrt(1+1/sv(i)^2);
    cs(i) = sn(i)/sv(i);
  end
end
end