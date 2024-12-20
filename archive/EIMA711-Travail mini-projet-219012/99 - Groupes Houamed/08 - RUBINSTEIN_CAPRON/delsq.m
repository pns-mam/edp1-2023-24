function D = delsq(G)
[m,n] = size(G);

% Indices of interior points
p = find(G);

% Connect interior points to themselves with 4's.
i = G(p);
j = G(p);
s = 4*ones(size(p));

% for k = north, east, south, west
for k = [-1 m 1 -m]
   % Possible neighbors in k-th direction
   Q = G(p+k);
   % Index of points with interior neighbors
   q = find(Q);
   % Connect interior points to neighbors with -1's.
   i = [i; G(p(q))];
   j = [j; Q(q)];
   s = [s; -ones(length(q),1)];
end
D = sparse(i,j,s);
