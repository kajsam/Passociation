function [Z, thresh] = cell_association_matrix(X, A, cell_effect, fig_nr)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation
%           cell_effect

n = size(X,1);

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    disp('Wrong dimensions')
    return
  end
else
  disp('Binary, please')
  return
end

Au = unique(A, 'rows', 'stable'); % Finding replicates. 
nu = size(Au,1);

% The association matrix
Z = zeros(n,nu);
for i = 1: nu
  Arow = repmat(Au(i,:),n,1);            % For each row, compared to all rows in X
  AndX = Arow & X;                      % Which entries are both 1
  sumA = sum(Au(i,:));
  sumA = max(sumA,1);                   % In case sumA = 0
  Z(:,i) = sum(AndX,2)./sumA;    % Sum up row-wise, normalised by the total number of 1's in the row 
end

Cell = repmat(cell_effect,1,nu);

figure(fig_nr), subplot(1,3,1), imagesc(Z), colormap(gray)
subplot(1,3,3), imagesc(Cell), colormap(gray)

Z = Z - Cell;
subplot(1,3,2), imagesc(Z), colormap(gray)

thresh = zeros(1,nu);
for i = 1: nu
  z = Z(:,i);
  thresh(i) = otsu_thresh(Z(:,i),256);
  z(z<thresh(i)) = 0;
  Z(:,i) = logical(z);
end

Z = logical(Z);
Z = unique(Z', 'rows', 'stable'); % Finding replicates. 
Z = Z';
figure(fig_nr+1), imagesc(Z), colormap(gray)

