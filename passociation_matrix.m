function Z = passociation_matrix(X,A,tau)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation

n = size(X,1);


if islogical(X) && islogical(A)
  if all(size(X)==size(A))
    'Welcome'
  else
    'Wrong dimensions'
    return
  end
else
  'Binary, please'
  return
end

Z = zeros(n,n);
for i = 1: n
  Arow = repmat(A(i,:),n,1);
  AndX = Arow & X;
  Z(:,i) = sum(AndX,2)./sum(A(i,:));
end

Z(Z<tau) = 0;

Z = logical(Z);

