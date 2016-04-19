X = [0.8;0.5];
Y = [1.3];
eta = 10

for i=1:999
	x = 5*(2*rand(2,1)-1);
	X = [X x];
	Y = [Y;sum(x)];
end
[m,n]=size(X);
W0 = zeros(m,1);

W1 = W0;
W2 = zeros(m,1);
Y1 = zeros(n,1);
num_iter = 1;
while (num_iter < 10000 & norm(Y1-Y) > 0.001)
	for i=1:n
		x = X(:,i);
		W2 = W1;
		x = X(:,i);
		z = dot(W2,x);
		Y1(i) = transfer(z,'lin');
		W1 = W2 + eta*(Y(i)-Y1(i))*x;
		W1 = W1/norm(W1);
		norm(Y-Y1);
		num_iter = num_iter+1;
		if norm(Y-Y1) < 0.001
			break
		endif
	end
	if norm(Y-Y1) < 0.001
		break
	endif
	if norm(W1-W2) < 0.01
		break
	endif
end
norm(Y1-Y)
W2
%abs(Y1-Y)

