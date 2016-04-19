x1 = [0.8;0.5;0];
x2 = [0.9;0.7;0.3];
x3 = [1.0;0.8;0.5];
x4 = [0;0.2;0.3];
x5 = [0.2;0.1;1.3];
x6 = [0.2;0.7;0.8];
X=[x1 x2 x3 x4 x5 x6];
eta = 1

[m,n]=size(X)
W0 = zeros(m,1);
Y = [-1;-1;-1;1;1;1]

W1 = W0
W2 = zeros(m,1)
Y1 = zeros(n,1)
while !(all(Y1==Y)==1)
	for i=1:n
		x = X(:,i)
		W2 = W1;
		x = X(:,i)
		z = dot(W2,x)
		Y1(i) = transfer(z,'bipolar')
		W1 = W2 + eta*(Y(i)-Y1(i))*x
		norm(W1-W2);
		if all(Y1==Y)==1
			break
		endif
	end
end
W2
norm(W1-W2)

