function Z = passociation_matrix(X,A,tau)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation
%           tau - threshold, scalar

n = size(X,1);

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    'Wrong dimensions'
    return
  end
else
  'Binary, please'
  return
end

% The association matrix
Z = zeros(n,n);
for i = 1: n
  Arow = repmat(A(i,:),n,1);            % For each row, compared to all rows in X
  AndX = Arow & X;                      % Which entries are both 1
  Z(:,i) = sum(AndX,2)./sum(A(i,:));    % Sum up row-wise, normalised by the total number of 1's in the row
end

Z(Z<tau) = 0;                           % Binarise according to tau

% Finding replicates. *************************
% ********************** Can possibly be done more elegantly
Z = logical(Z);
figure(10), subplot(1,2,1), imagesc(Z)
Zdiff = false(n,n);
for i = 1: n
  Zcol = repmat(Z(:,i),1,n);
  Zxor = xor(Zcol,Z);
  Zdiff(i,:) = logical(sum(Zxor,1));
end

Zeq = ~Zdiff;
UpTr = triu(Zeq,1);
Zeq = Zeq & UpTr;
[~,col] = find(Zeq);

Z(:,col) = [];                          % Delete replicates

figure(1), imagesc(Z)

