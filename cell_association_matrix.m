function [binZ, thresh] = cell_association_matrix(X, A, cell_effect, fig_nr)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation
%           cell_effect

[n,d] = size(X);

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
sum0 = logical(sum(Au,2));                   % Column of all zeros
Au = Au(sum0,:);
nu = size(Au,1);

% The association matrix
Z = zeros(n,nu);
XnA = zeros(1,n);
for i = 1: nu
  for l = 1 :n
    XnA(l) = sum(X(l,:) & Au(i,:));
  end
  sumA = sum(Au(i,:));
  sumA = max(sumA,1);                   % In case sumA = 0
  Z(:,i) = XnA./sumA;    % Sum up row-wise, normalised by the total number of 1's in the row 
end

Cell = repmat(cell_effect,1,nu);

if length(fig_nr) > 1 
  figure(fig_nr(2)), subplot(1,3,1), imagesc(Z), colormap(gray), title('Association matrix')
  subplot(1,3,3), imagesc(Cell), colormap(gray), title('Cell effect')
end
Z = Z - Cell;

num_bins = 256; 
thresh = zeros(1,nu);
binZ = false(nu,n);
for i = 1: nu
  z = Z(i,:);
  thresh(i) = otsu_thresh(Z(i,:),num_bins);
  z(z<thresh(i)) = 0;
  binZ(i,:) = logical(z);
end

binZ = unique(binZ', 'rows', 'stable'); % Finding replicates. 
binZ = binZ';

figure(fig_nr(1)), imagesc(binZ), colormap(gray), 
title(strcat(num2str(size(Z,2)), ' candidate columns'))

if length(fig_nr) > 1
  figure(fig_nr(2)),subplot(1,3,2), imagesc(Z), colormap(gray), title('Cell-adjusted association')
  subplot(1,3,1)
end