x1 = [1;-2.0;1.5;0.0];
x2 = [1.0;-0.5;-2.0;-1.5];
x3 = [0.0;1.0;-1.0;1.5];
X=[x1 x2 x3];
eta = 0.2

W0 = [1.0;-1.0;0.0;0.5];
Y = [-1;-1;1]

W1 = W0
W2 = zeros(4,1)
Y1 = zeros(3,1)
while !(all(Y1==Y)==1)
	for i=1:3
		x = X(:,i)
		W2 = W1;
		x = X(:,i)
		z = dot(W1,x)
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

