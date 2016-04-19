x1 = [0.8;0.5;0];
x2 = [0.9;0.7;0.3];
x3 = [1.0;0.8;0.5];
x4 = [0;0.2;0.3];
x5 = [0.2;0.1;1.3];
x6 = [0.2;0.7;0.8];
X=[x1 x2 x3 x4 x5 x6];
eta = 50

[m,n]=size(X);
W0 = zeros(m,1);
Y = [-1;-1;-1;1;1;1];

W1 = W0;
W2 = zeros(m,1);
Y1 = zeros(n,1);
num_iter = 1;
while (num_iter < 500 & norm(Y1-Y) > 0.01)
	for i=1:n
		x = X(:,i);
		W2 = W1;
		x = X(:,i);
		z = dot(W2,x);
		Y1(i) = 2*transfer(z,'sigmoid')-1
		W1 = W2 + eta*(Y(i)-Y1(i))*x*0.5*(1-Y1(i)*Y1(i))
		norm(Y-Y1)
		num_iter = num_iter+1;
		if norm(Y-Y1) < 0.01
			break
		endif
	end
	if norm(Y-Y1) < 0.01
		break
	endif
end
Y1
norm(Y1-Y)

