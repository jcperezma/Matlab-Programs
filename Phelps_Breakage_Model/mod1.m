function answer = mod1(X,Y)
  tol = 1e7;
  X = round(X*tol)/tol;
  Y = round(Y*tol)/tol;
  answer = builtin('mod',X,Y);
end