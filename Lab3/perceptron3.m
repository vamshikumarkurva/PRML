x1 = [1;1];
x2 = [1;-1];
x3 = [-1;1];
x4 = [-1;-1];
X=[x1 x2 x3 x4];
eta = 1

[m,n]=size(X)
W0 = zeros(m,1);
Y = [1;-1;-1;-1]

W1 = W0
W2 = zeros(m,1)
Y1 = zeros(n,1)
num_iter = 1
while !(all(Y1==Y)==1)
	for i=1:n
		x = X(:,i)
		W2 = W1
		x = X(:,i);
		z = dot(W2,x)
		Y1(i) = transfer(z,'bipolar')
		W1 = W2 + eta*(Y(i)-Y1(i))*x
		W2 = W2/norm(W2)
		norm(W1-W2)
		if all(Y1==Y)==1
			break
		endif
	end
	num_iter = num_iter+1
	if num_iter > 100
		break
	endif
	if all(Y1==Y)==1
		break
	endif
end
W2
norm(W1-W2)

