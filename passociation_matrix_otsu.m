function [Z, thresh] = passociation_matrix_otsu(X, A, smooth)

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
  sumA = sum(A(i,:));
  sumA = max(sumA,1);                   % In case sumA = 0
  Z(:,i) = sum(AndX,2)./sumA;    % Sum up row-wise, normalised by the total number of 1's in the row 
end

if smooth
  smoothZ = Z;
  smoothZ(logical(eye(size(Z)))) = NaN;
  Z = smoothdata(smoothZ,2,'movmedian',3);
  Z(logical(eye(size(Z)))) = 1;
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

% figure, imagesc(Z), colormap(gray), title(median(thresh))

% Finding replicates. *************************
% ********************** Can possibly be done more elegantly
% *** Faster for sure


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

Z(:,unique(col)) = [];                          % Delete replicates

