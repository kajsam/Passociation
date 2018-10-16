function [Z, thresh] = thresh_matrix(Z, fig_nr)

% Kajsa Mollersen (kajsa.mollersen@uit.no) October 15th 2018

% Input:    X (n x d)           - binary gene expression matrix
%           A (n x d)           - binary structure matrix approximation
%           alphabeta           - the Bernoulli pi parameters

% This function calculates n column vectors of length n, that will be
% candidates column vectors for the next rank K approximation A of X.

size(Z)
n = size(Z,1);

% if fig_nr
%   figure(fig_nr), subplot(1,2,1), imagesc(Z), colormap(gray)
% end

num_bin = 256
thresh = zeros(1,n);
step = (max(Z(:)) - min(Z(:)))/100;
Zt = false(size(Z)); 
for t = min(Z(:))+step:step:max(Z(:))-step
  for i = 1: n
    z = Z(:,i);
%   z = z - min(z);
%   z = z./max(z);
%   thresh(i) = otsu_thresh(Z(:,i),num_bin);
%   figure(fig_nr+1), subplot(1,2,1), histogram(Z(:,i),num_bin)
%   subplot(1,2,2), histogram(z,num_bin), title(thresh(i))
    z(z<t) = 0;
   
    Zt(:,i) = logical(z); 
%   pause
  end
  
  figure(fig_nr), imagesc(Zt), colormap(gray), title(t)
  pause
end

if fig_nr  
  figure(fig_nr), subplot(1,2,2), imagesc(Z), colormap(gray), title(median(thresh))
end
% Finding replicates. *************************
% ********************** Can possibly be done more elegantly
% *** Faster for sure

Zdiff = false(size(Z));
for i = 1: n
  Zcol = repmat(Z(:,i),1,size(Z,2));
  Zxor = xor(Zcol,Z);
  Zdiff(i,:) = logical(sum(Zxor,1));
end

Zeq = ~Zdiff;
UpTr = triu(Zeq,1);
Zeq = Zeq & UpTr;
[~,col] = find(Zeq);

Z(:,unique(col)) = [];        
% Delete replicates
