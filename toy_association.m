function toy_association

Y = [1 1 0 0;
     1 1 0 1;
     0 1 0 0];
 
Y = logical(Y);

tau = 0.9;

Z = passosiation_matrix(Y,Y,0.9)

[W, H] = asso(double(Y'), 3, 0.9);

W'