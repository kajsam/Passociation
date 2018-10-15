function [Z, thresh, em] = prociation_matrix(X, A, cell_effect)

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

Cell = repmat(cell_effect,1,n);

figure, subplot(1,2,1), imagesc(Z(1:100,:)), colormap(gray)
subplot(1,2,2), imagesc(Cell(1:100,:)), colormap(gray)

sum0 = sum(Z,1);
Z(:,~sum0) = [];
figure, subplot(1,2,1), imagesc(Z, [0 1]), colormap(gray)
Cell(:,~sum0) = [];
Z = Z - Cell;
subplot(1,2,2), imagesc(Z, [0 1]), colormap(gray)
tau = otsu_thresh(Z,256);
Z0 = Z;
step = (max(Z(:)) - min(Z(:)))/20
for t = min(Z(:))+step:step: max(Z(:))-step
    Z = Z0;
  Z(Z<t) = 0;
Z = logical(Z);
sum0 = sum(Z,1);
Z(:,~sum0) = [];
figure, imagesc(Z), colormap(gray), title(t)
pause
end
  
Z = Z0;

num_bin = floor(n/10);
thresh = zeros(1,size(Z,2));
em = zeros(1,size(Z,2));
for i = 1: size(Z,2)
  z = Z(:,i);
  [thresh(i) em(i)] = otsu_thresh(Z([1:i-1 i+1:size(Z,2)],i),num_bin);
  z(z<thresh(i)) = 0;
  % Z(:,i) = logical(z);
end
Z(Z<tau) = 0;
Z = logical(Z);
sum0 = sum(Z,1);
Z(:,~sum0) = [];
figure, imagesc(Z), colormap(gray), title(median(thresh))

% Finding replicates. *************************
% ********************** Can possibly be done more elegantly
% *** Faster for sure

Zdiff = false(size(Z));
for i = 1: size(Z,2)
  Zcol = repmat(Z(:,i),1,size(Z,2));
  Zxor = xor(Zcol,Z);
  Zdiff(i,:) = logical(sum(Zxor,1));
end

Zeq = ~Zdiff;
UpTr = triu(Zeq,1);
Zeq = Zeq & UpTr;
[~,col] = find(Zeq);

Z(:,unique(col)) = [];                          % Delete replicates
figure, imagesc(Z), colormap(gray)
