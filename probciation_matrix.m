function Z = probciation_matrix(X, A, alphabeta, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) October 15th 2018

% Input:    X (n x d)           - binary gene expression matrix
%           A (n x d)           - binary structure matrix approximation
%           alphabeta           - the Bernoulli pi parameters

% This function calculates n column vectors of length n, that will be
% candidates column vectors for the next rank K approximation A of X.

n = size(X,1);

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    disp('Wrong dimensions')
    return
  end
else
  disp('Logical, please')
  return
end

alpha = alphabeta(1,:);
beta = alphabeta(2,:);

% The probability/likelihood matrix

% This is obviously not what I want!
Z = zeros(n,n);
for i = 1: n
  for ii = 1: n  % Compare column i in X to column ii in Pi
    n00 = sum(~X(i,:) & ~A(ii,:));
    n01 =sum(~X(i,:) & A(ii,:));
    n10 = sum(X(i,:) & ~A(ii,:));
    n11 = sum(X(i,:) & A(ii,:));
    
    Z(ii,i) = (n10+n11)*alpha(ii) + n11*beta(ii) ...
        - (n00+n10)*log(1 + exp(alpha(ii))) ...
        - (n01 + n11)*log(1 + exp(alpha(ii) + beta(ii)));
  end
end