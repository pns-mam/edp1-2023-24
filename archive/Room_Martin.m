function U=RoomTemperature(ot,dt,ht,n);
% ROOMTEMPERATURE computes the room temperature in our living room
%   U=RoomTemperature(ot,dt,ht,n); takes the outside temperature ot,
%   the door temperature dt and the heater temperature ht and the
%   number of gridpoints n and then computes the room temperature in
%   our living room

x=linspace(-1,1,n);                              % generate grid
[X,Y]=meshgrid(x,x);
G=((X>-1) & (X<0.6) & (Y>-0.5) & (Y<0.8)) | ...  % area of our living room
  ((X>-0.2) & (X<0.6) & (Y>-1) & (Y<-0.5)) | ... 
  ((X>-0.6) & (X<0) & (Y>0.8) & (Y<1));
H=((X>-0.6) & (X<0) & (Y>0.5) & (Y<0.75));       % original heater location
%H=((X>-0.6) & (X<0) & (Y>-0.1) & (Y<0.15));     % heater in the center
%H=((X>0.4) & (X<0.6) & (Y>-0.3) & (Y<0.3));     % heater on the wall
k=find(G);
G=zeros(size(G));                   % Convert from logical to double
G(k)=(1:length(k))';                % Indices for the matrix
A=delsq(G);                         % Laplacian in the interior
do=[];                              % door indices
for i=2:n-1,                        % add Neumann conditions for insulated
  for j=2:n-1,                      % walls
    no=G(i,j);
    if no~=0,
      if G(i,j-1)==0,               % Neumann condition on the left wall
        A(no,no)=A(no,no)-1;
      end;
      if G(i,j+1)==0,               % Neumann condition on the right wall
        A(no,no)=A(no,no)-1;
      end;
      if G(i+1,j)==0 & i<n-1,       % keep Dirichlet conditions for window 
        A(no,no)=A(no,no)-1;
      end;
      if G(i-1,j)==0,
        if (X(i,j)>-0.8 & X(i,j)<-0.3), % keep Dirichlet conditions for door
          do=[do no];
        else
          A(no,no)=A(no,no)-1;
        end;
      end;  
    end;
  end;
end;
h=2/(n-1);
A=-A/h^2;                           % scale the Laplacian with h
b=zeros(length(k),1);
wi=G(end-1,find(G(end-1,:)>0));     % find window indices
he=G(find(H));                      % find heater indices
b(wi)=-1/h^2*ot;                    % add heating and Dirichlet conditions
b(he)=-ht;
b(do)=-1/h^2*dt;
u=A\b;                              % solve using sparse reordered LU
U=G;
U(G>0)=full(u(G(G>0)));             % put solution onto the grid
mesh(X,Y,U);
axis('ij');
