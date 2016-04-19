x1 = [1;-2.0;1.5;0.0];
x2 = [1.0;-0.5;-2.0;-1.5];
x3 = [0.0;1.0;-1.0;1.5];
W0 = [1.0;-1.0;0.0;0.5];
X=[x1 x2 x3];
eta = 0.2
W1 = W0/norm(W0);
W2 = zeros(4,1);

num_iter = 1
while norm(W1-W2) > 0.01
	for i=1:3
		W2 = W1;
		x = X(:,i);
		z = dot(W1,x);
		y = 2*transfer(z,'sigmoid')-1
		W1 = W2 + eta*y*x;
		if norm(W1)~=0
			W1 = W1/norm(W1);
		endif
		num_iter = num_iter+1;
		norm(W1-W2);
		if num_iter > 60
			break
		endif
		if norm(W1-W2)<0.01
			break
		endif
	end
	#{
	if num_iter > 500
		break
	endif
	#}
end
W2
norm(W1-W2)
