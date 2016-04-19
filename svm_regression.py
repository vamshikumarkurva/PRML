import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
import math


GAUSSIAN_KERNEL = 1
LINEAR_KERNEL = 0
POLYNOMIAL_KERNEL = 0
def kernel(x,y):
	if GAUSSIAN_KERNEL:
		z = x-y
		gamma = 120 # gamma = 1/(2*sigma*sigma)
		K = np.exp(-gamma*np.power(LA.norm(z),2))
		return K
	elif LINEAR_KERNEL:
		K=x*y
		return K.sum()
	elif POLYNOMIAL_KERNEL:
		K=x*y
		k1 = K.sum()
		order = 1
		k = math.pow(k1+2,order)
		return k

def SVMoutput(data,alpha,point,b):
	m,n = data.shape
	out = 0
	for i in range(m):
		out += alpha[i]*kernel(data[i,0:n-1],point)
	out += b[0]
	return out

def findError(alpha,data,b,C,output,EPS):
	m,n=data.shape
	for i in range(m):
		y2 = data[i][n-1]
		point = data[i,0:n-1]
		output[i] = SVMoutput(data,alpha,point,b)
		error[i] = output[i]-y2
		output[i] = abs(error[i])
	return error 

def DualobjFunc(alpha,data,EPS):
	out1 = 0
	out2 = 0
	m,n=data.shape
	for i in range(m):
		for j in range(m):
			y1 = data[i][n-1]
			y2 = data[j][n-1]
			point1 = data[i,0:n-1]
			point2 = data[j,0:n-1]
			out1 += alpha[i]*alpha[j]*kernel(point1,point2)
	Y = data[:,n-1]
	out2 = alpha*Y
	out = out2.sum()-EPS*(abs(alpha).sum())-0.5*out1
	return out

def PrimalobjFunc(alpha,data,EPS,b,C):
	out1 = 0
	out2 = 0
	m,n=data.shape
	for i in range(m):
		for j in range(m):
			y1 = data[i][n-1]
			y2 = data[j][n-1]
			point1 = data[i,0:n-1]
			point2 = data[j,0:n-1]
			out1 += alpha[i]*alpha[j]*kernel(point1,point2)
	E1 = np.zeros(m)
	E2 = np.zeros(m)
	for i in range(m):
		y2 = data[i][n-1]
		point = data[i,0:n-1]
		out = SVMoutput(data,alpha,point,b)
		E1[i] = max(0,out-y2-EPS)
		E2[i] = max(0,y2-out-EPS)
	out = 0.5*out1+C*(E1.sum()+E2.sum())
	return out

def takeStep(i1,i2,alpha,data,b,C,error,output,EPS):
	eps = 0.001
	if i1==i2:
		return 0
	alph1 = alpha[i1]
	alph2 = alpha[i2]
	y1 = data[i1][n-1]
	y2 = data[i2][n-1]
	point1 = data[i1,0:n-1]
	point2 = data[i2,0:n-1]
	error = findError(alpha,data,b,C,output,EPS)
	E1 = error[i1]
	E2 = error[i2]
	if alph2 > 0:
		L_2 = 0
		H_2 = C
	else:
		L_2 = -C
		H_2 = 0
	if alph1 > 0:
		L_1 = 0
		H_1 = C
	else:
		L_1 = -C
		H_1 = 0
	L = max(L_2, alph1+alph2-H_1)
	H = min(H_2, alph1+alph2-L_1)
	if alph1==0 and alph2==0:
		L = max(-C, alph1+alph2-C)
		H = min(C, alph1+alph2+C)
	if L==H:
		return 0
	k11 = kernel(point1, point1)
	k12 = kernel(point1, point2)
	k22 = kernel(point2, point2)
	eta = 2*k12-k11-k22
	if eta < 0:
		a2 = alph2 - (E1-E2-EPS*(np.sign(alph2)-np.sign(alph1)))/eta
		if a2 < L:
			a2 = L
		elif a2 > H:
			a2 = H
	else:
		alpha1 = alpha
		alpha1[i2] = L
		Lobj = DualobjFunc(alpha1,data,EPS)
		alpha1[i2]=H
		Hobj = DualobjFunc(alpha1,data,EPS)
		if Lobj > Hobj+eps:
			a2 = L
		elif Lobj < Hobj-eps:
			a2 = H
		else:
			a2 = alph2
	if abs(a2-alph2) < eps*(a2+alph2+eps):
		return 0
	a1 = alph1+alph2-a2
	alpha[i1]=a1
	alpha[i2]=a2

	b1 = b[0]-E1+k11*(a1-alph1)+k12*(a2-alph2)
	b2 = b[0]-E2+k12*(a1-alph1)+k22*(a2-alph2)
	# updating threshold
	if b1==b2:
		b[0]=b1
	else:
		if abs(a1) < C and abs(a1)!=0:
			b[0] = b1
		elif abs(a2) < C and abs(a2)!=0:
			b[0] = b2
		else:
			b[0] = (b1+b2)/2
	# end
	return 1

def secondHueristic(alpha,data,b,c,error,i2,output,EPS):
	m,n=data.shape
	error = findError(alpha,data,b,C,output,EPS)
	E2 = error[i2]
	maxi = None
	max_ind = None
	for i in range(m):
		if maxi is None:
			maxi = abs(error[i]-E2)
			max_ind = i
		else:
			if abs(error[i]-E2) > maxi:
				maxi = abs(error[i]-E2)
				max_ind = i
	return max_ind 
	

def examineExample(i2,alpha,data,b,C,error,output,EPS,primal,dual):
	tol = 0.001
	m,n = data.shape
	error = findError(alpha,data,b,C,output,EPS)
	alph2 = alpha[i2]
	y2 = data[i2][n-1]
	E2 = error[i2]
	x = DualobjFunc(alpha,data,EPS)
	y = PrimalobjFunc(alpha,data,EPS,b,C)
	print x,y
	dual.append(x)
	primal.append(y)
	if (E2 > -EPS+tol and alph2 > 0) or (E2 < EPS-tol and alph2 < 0 ) or (E2 > tol and alph2==0) or (E2 < -tol and alph2==0):
		# finding no of non-zero, non-C alphas
		ind = np.where((alpha > -C) & (alpha < C))[0]
		indlist = ind.tolist()
		k = len(indlist)
		if k!=0:
			i1 = secondHueristic(alpha,data,b,C,error,i2,output,EPS)
			if takeStep(i1,i2,alpha,data,b,C,error,output,EPS):
				return 1
		
		if k!=0:
			for i1 in indlist:
				if takeStep(i1,i2,alpha,data,b,C,error,output,EPS):
					return 1
		
		for i1 in range(m):
			if i1!=i2:
				if takeStep(i1,i2,alpha,data,b,C,error,output,EPS):
					return 1
	return 0

def plotdata(data):
	m,n = data.shape
	data_x=[]
	data_y=[]
	plt.figure(1)
	for i in range(m):
		data_x.append(data[i][0])
		data_y.append(data[i][1])
	plt.plot(data_x,data_y,'ro')
	plt.xlabel('------> x')
	plt.ylabel('------> y')
	plt.show()

if __name__ == '__main__':
	data = np.loadtxt("/home/iistlab/Desktop/vamshi/assign2/data3.csv",delimiter=",")
	(m,n)=data.shape
	plotdata(data)
	alpha = np.zeros(m)
	error = np.zeros(m)
	output = np.zeros(m)
	b = np.zeros(1)
	C = 0.2
	EPS = 0.001
	max_iter = 50
	num_iter=0

	numChanged = 0
	examineAll = 1
	primal=[]
	dual=[]
	while numChanged > 0 or examineAll:
		numChanged = 0
		if examineAll:
			for i in range(m):
				numChanged += examineExample(i,alpha,data,b,C,error,output,EPS,primal,dual)
		else:
			for i in range(m):
				if abs(alpha[i]) > 0 and abs(alpha[i]) < C:
					numChanged += examineExample(i,alpha,data,b,C,error,output,EPS,primal,dual)
		if examineAll == 1:
			examineAll = 0
		elif numChanged == 0:
			examineAll = 1
		print output
		ind = np.where(output<0.8)[0]
		num_iter += 1
		indlist = ind.tolist()
		k = len(indlist)
		print k
		l = 0.90*m
		print "iteration: ", num_iter
		if k > l or num_iter >= max_iter:
			break

	print "threshold: ", b
	print "Alpha: ", alpha
	print "output: ", output

	'''
	# alpha obtained for data4 
	b[0] = np.array([ 5.85697632])
	alpha=np.array([-0.2        -0.2         0.2        -0.2        -0.2         0.          0.2
	 -0.2         0.2         0.2         0.2        -0.2        -0.2        -0.2
	 -0.2        -0.2        -0.2         0.2         0.2        -0.2         0.2
	  0.2        -0.2        -0.2        -0.2         0.2         0.2         0.2
	  0.2         0.2         0.2        -0.2        -0.2         0.2        -0.2
	  0.2         0.2         0.2        -0.2        -0.2        -0.2        -0.2
	 -0.2         0.2        -0.2         0.2         0.2        -0.2        -0.2
	  0.2        -0.2         0.2        -0.2         0.2        -0.2        -0.2
	 -0.2         0.2         0.2         0.2         0.2        -0.2        -0.2
	 -0.2         0.2         0.2        -0.2         0.2        -0.2         0.2
	 -0.2        -0.2         0.2         0.2        -0.2        -0.2         0.2
	 -0.2         0.2        -0.2         0.2        -0.2         0.2        -0.2
	  0.2         0.2        -0.2         0.2        -0.2         0.2        -0.2
	 -0.2        -0.2         0.2        -0.2         0.2        -0.2        -0.2
	 -0.2         0.2        -0.2        -0.2        -0.2        -0.18928583
	  0.2        -0.2        -0.2         0.2         0.2        -0.2         0.2
	  0.2        -0.2         0.2        -0.2         0.2         0.2        -0.2
	 -0.2        -0.2         0.2         0.2        -0.2        -0.2         0.2
	  0.2        -0.2         0.2        -0.2        -0.2         0.2        -0.2
	 -0.2         0.2        -0.2         0.2        -0.2         0.2         0.2
	 -0.2         0.19366612  0.2         0.2        -0.2         0.2        -0.2
	  0.2        -0.2         0.2        -0.2         0.2         0.2        -0.2
	 -0.15205722 -0.2         0.2        -0.2         0.14397241 -0.2         0.2
	 -0.2         0.2         0.04531343  0.2        -0.2         0.2         0.2
	  0.2        -0.2        -0.2        -0.08305153  0.2         0.2         0.2
	 -0.2         0.2         0.2        -0.2        -0.2         0.2         0.2
	  0.2        -0.2        -0.2         0.2        -0.14486836  0.2        -0.052522
	  0.2         0.2         0.2        -0.16116701  0.2         0.2        -0.2
	 -0.2         0.2         0.2        -0.2         0.2       ])
	'''
	plt.figure(2)
	for i in range(m):
		if alpha[i]==0:
			plt.plot(data[i][0],data[i][1],'ro')
		else:
			plt.plot(data[i][0],data[i][1],'bo')

	# plotting the function
	xmin = min(data[:,0])
	xmax = max(data[:,0])
	boundary_x = []
	boundary_y = []
	for i in np.linspace(xmin-1, xmax+1, 400):
		point = np.array([i])
		out = SVMoutput(data,alpha,point,b)
		boundary_x.append(i)
		boundary_y.append(out)
	plt.plot(boundary_x,boundary_y,'r--') 
	plt.xlabel('------> x')
	plt.ylabel('------> y')
	plt.show()

	k = min(len(primal), len(dual))
	plt.figure(3)
	plt.plot(range(k),primal[:k],'r-',label='Primal obj func')
	plt.plot(range(k),dual[:k],'b-',label='Dual obj func')
	plt.xlabel('------> iterations')
	plt.ylabel('------> objective function')
	plt.title('convergence plot')
	plt.axis([0,k,min(min(primal),min(dual)),max(max(primal),max(dual))])
	plt.show()
	
	
		
