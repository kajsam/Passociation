function Zredux = reduxZoverlap(Z, min_class)

if ~islogical(Z)
  disp('Logical, please')
  return
end

[n, nZ] = size(Z);
sumZ = sum(Z,1);
[~, idx] = sort(sumZ,'descend');
Z = Z(:,idx);
   
keep = false(1,nZ);
for k = 1 : nZ
  Zcand = Z(:,k+1:end);
  del = false(1,nZ);
  for i = k+1:nZ
    if sum(Z(~Z(:,k),i)) > min_class
      del(i) = true;
    end
  end
  del(1:k) = [];
  Zcand(:,del) = [];
   
  repr = logical(sum(Zcand,2));
  
  if sum(Z(:,k) & repr) < sum(Z(:,k))
    keep(k) = true;
  end
end
Zredux = Z(:,keep);