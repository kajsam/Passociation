function [Z, thresh] = passociation_matrix_otsu(X, A, smooth, fig_nr)

% Input:    X - binary gene expression matrix
%           A - binary structure matrix approximation
%           smooth - option for smoothin, not a good idea after all

[n,d] = size(X);

if islogical(X) && islogical(A)
  if ~all(size(X)==size(A))
    disp('Wrong dimensions')
    return
  end
else
  disp('Logical, please')
  return
end

% The association matrix
Z = zeros(n,n);
AndX = zeros(n,d);
for i = 1: n
  for j = 1 :n
    AndX(j,:) = X(j,:)& A(i,:);
  end
  % Arow = repmat(A(i,:),n,1);            % For each row, compared to all rows in X
  % AndX = Arow & X;                      % Which entries are both 1
  sumA = sum(A(i,:));
  sumA = max(sumA,1);                   % In case sumA = 0
  Z(:,i) = sum(AndX,2)./sumA;    % Sum up row-wise, normalised by the total number of 1's in the row 
end

if smooth % Do not use this
  smoothZ = Z;
  smoothZ(logical(eye(size(Z)))) = NaN;
  Z = smoothdata(smoothZ,2,'movmedian',3);
  Z(logical(eye(size(Z)))) = 1;
end

if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
end 

num_bin = floor(n/10);
thresh = zeros(1,n);
for i = 1: n
  z = Z(:,i);
  thresh(i) = otsu_thresh(Z([1:i-1 i+1:n],i),num_bin);
  z(z<thresh(i)) = 0;
  Z(:,i) = logical(z);
end

Z = logical(Z);

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(Z), colormap(gray), title(median(thresh))
end
Z = unique(Z', 'rows', 'stable'); % Finding replicates. 
Z = Z';

