function [binZ, thresh] = association_matrix_cell_adjusted(X, A, pi1, min_class, fig_nr)

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
  disp('Logical, please')
  return
end

Au = unique(A, 'rows', 'stable'); % Finding replicates. 
sum0 = logical(sum(Au,2));                   % Column of all zeros
Au = Au(sum0,:);
nu = size(Au,1)

% The association matrix
Z = zeros(n,nu);
XnA = zeros(1,n);
nXA = zeros(1,n);
logpi1 = log(pi1);
logpi1compl = log(ones(1,n)- pi1);
for i = 1: nu
  for l = 1 :n
    XnA(l) = sum(X(l,:) & Au(i,:));
    nXA(l) = sum(~X(l,:) & Au(i,:));
  end
  
  Z(:,i) = XnA.*logpi1 + nXA.*logpi1compl;    % Sum up row-wise, normalised by the total number of 1's in the row 
i
end

if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
  title('Cell effect adjusted association matrix')
end 

num_bins = 256; 
thresh = zeros(1,nu);
binZ = false(n,nu);
for i = 1: n
  z = Z(i,:) - min(Z(i,:));
  z = z./max(z);
  % figure(fig_nr+1), hist(z,num_bins)
  thresh(i) = otsu_thresh(z,num_bins);
  % title(thresh(i))
  % pause
  z(z<thresh(i)) = 0;
  binZ(i,:) = logical(z);
end

binZ = unique(binZ', 'rows', 'stable'); % Finding replicates. 
binZ = binZ';

sumZ = sum(binZ,1);

idx0 = sumZ < min_class;    % Smaller than the smallest class size
binZ(:,idx0) = [];          % Noise, delete it

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(binZ) 
  colormap(gray),title(strcat(num2str(size(binZ,2)), ' candidate columns'))
end
