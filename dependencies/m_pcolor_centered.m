% PCOLOR_CENTERED(first 3 basic pcolor arguments)
% 
% Makes pcolor actually plot the matrix you give it (e.g. give it a 10x10
% matrix, plot a 10x10 matrix, not a 9x9 matrix).

function PC = m_pcolor_centered(X,Y,D)

dx = X(1,2) - X(1,1);
X_ = [X , X(:,end) + dx ; ...
      X(end,:) , X(end,end) + dx];

dy = Y(2,1) - Y(1,1);
Y_ = [Y , Y(:,end) ; ...
      Y(end,:) + dy , Y(end,end) + dy];

D_ = [D , nan(size(D,1),1) ; ...
      nan(1,size(D,2)) , nan ];

PC = m_pcolor(X_ - dx/2, Y_ - dy/2, D_);

end
