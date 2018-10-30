function [binZ, thresh] = association_matrix_gene_adjusted(X, A, pi1, min_class, fig_nr)

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
sum0 = sum(Au,2);
Au = Au(logical(sum0),:);
nu = size(Au,1);

% The association matrix
Z = zeros(n,nu);
XnA = zeros(n,nu);
nXA = zeros(n,nu);
logpi1 = log(pi1);
logpi1compl = log(ones(d,1)- pi1);
for i = 1: nu
  for l = 1 :n
    XnA(l,i) = sum(logpi1(X(l,:) & Au(i,:)));
    nXA(l,i) = sum(logpi1compl(~X(l,:) & Au(i,:)));
  end
end
denom = sum(XnA + nXA,2);
for i = 1: nu
  Z(:,i) = (XnA(:,i) + nXA(:,i))./denom;    % Sum up row-wise
end

if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
  title('Gene effect adjusted association matrix')
end 
num_bins = 256; 
thresh = zeros(1,nu);
binZ = false(n,nu);
for i = 1: nu
  z = Z(:,i) - min(Z(:,i));
  z = z./max(z);
  % figure(fig_nr+1), hist(z,num_bins)
  thresh(i) = otsu_thresh(z,num_bins);
  % title(thresh(i))
  % pause
  z(z>thresh(i)) = 0;
  binZ(:,i) = logical(z);
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

