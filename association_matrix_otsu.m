function [binZ, thresh] = association_matrix_otsu(X, fig_nr)

% Input:    X - binary gene expression matrix

[n,d] = size(X);

if ~islogical(X)
  disp('Logical, please')
  return
end

% The association matrix
Z = zeros(n,n);
XnX = zeros(1,n);
for i = 1: n                        % For each row, compare to all other rows 
  for l = 1 :n
    XnX(l) = sum(X(l,:)& X(i,:));      % Which entries are both 1
  end           
  sumX = sum(X(i,:));
  sumX = max(sumX,1);               % In case sum = 0
  Z(:,i) = XnX./sumX;        % Sum up row-wise, normalised by the total number of 1's in the row 
end

Z(logical(eye(n))) = median(Z(:));
if fig_nr
  figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
  title('Association matrix')
end 

% Calculate the threshold for binarization. The diagonal is excluded.  
num_bins = 256; %floor(n/25);
thresh = zeros(1,n);
binZ = false(n,n);
for i = 1: n
  z = Z(i,:);
  % figure(fig_nr+1), hist(z([1:i-1 i+1:n]),num_bin)
  thresh(i) = otsu_thresh(Z(i,[1:i-1 i+1:n]),num_bins);
  % title(thresh(i))
  % pause
  z(z<thresh(i)) = 0;
  binZ(i,:) = logical(z);
end

binZ = unique(binZ', 'rows', 'stable'); % Finding replicates. 
binZ = binZ';

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(binZ) 
  colormap(gray), title(median(thresh))
end


