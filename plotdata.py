import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

GAUSSIAN_KERNEL = 1
LINEAR_KERNEL = 0
POLYNOMIAL_KERNEL = 0
def kernel(x,y):
	if GAUSSIAN_KERNEL:
		z = x-y
		gamma = 100 #gamma = 1/(2*sigma*sigma)
		K = np.exp(-gamma*np.power(LA.norm(z),2))
		return K
	elif LINEAR_KERNEL:
		K=x*y
		return K.sum()
	elif POLYNOMIAL_KERNEL:
		K=x*y
		order = 3
		k = np.power(K.sum(),order)
		return k

def SVMoutput(data,alpha,point,b):
	m,n = data.shape
	out = 0
	for i in range(m):
		out += alpha[i]*data[i][n-1]*kernel(data[i,0:n-1],point)
	out += b[0]
	return out

def findError(alpha,data,b,C,output):
	m,n=data.shape
	for i in range(m):
		y2 = data[i][n-1]
		point = data[i,0:n-1]
		error[i] = SVMoutput(data,alpha,point,b)-y2
		output[i] = np.sign(SVMoutput(data,alpha,point,b)*y2)
	return error 

def DualobjFunc(alpha,data):
	out = 0
	m,n=data.shape
	for i in range(m):
		for j in range(m):
			y1 = data[i][n-1]
			y2 = data[j][n-1]
			point1 = data[i,0:n-1]
			point2 = data[j,0:n-1]
			out += alpha[i]*alpha[j]*y1*y2*kernel(point1,point2)
	out = alpha.sum()-0.5*out
	return out

def PrimalobjFunc(alpha,data,b,C,output):
	out = 0
	m,n=data.shape
	for i in range(m):
		for j in range(m):
			y1 = data[i][n-1]
			y2 = data[j][n-1]
			point1 = data[i,0:n-1]
			point2 = data[j,0:n-1]
			out += alpha[i]*alpha[j]*y1*y2*kernel(point1,point2)
	E = np.zeros(m)
	for i in range(m):
		y2 = data[i][n-1]
		point = data[i,0:n-1]
		E[i] = max(0,1-SVMoutput(data,alpha,point,b)*y2)
	out = 0.5*out + C*(E.sum())
	return out


def takeStep(i1,i2,alpha,data,b,C,error,output):
	eps = 0.001
	if i1==i2:
		return 0
	alph1 = alpha[i1]
	alph2 = alpha[i2]
	y1 = data[i1][n-1]
	y2 = data[i2][n-1]
	point1 = data[i1,0:n-1]
	point2 = data[i2,0:n-1]
	error = findError(alpha,data,b,C,output)
	E1 = error[i1]
	E2 = error[i2]
	s = y1*y2
	if s==1:
		L = max(0, alph1+alph2-C)
		H = min(C, alph1+alph2)
	else:
		L = max(0, alph1-alph2)
		H = min(C, C-alph1+alph2)
	if L==H:
		return 0
	k11 = kernel(point1, point1)
	k12 = kernel(point1, point2)
	k22 = kernel(point2, point2)
	eta = 2*k12-k11-k22
	if eta < 0:
		a2 = alph2 - y2*(E1-E2)/eta
		if a2 < L:
			a2 = L
		elif a2 > H:
			a2 = H
	else:
		alpha1 = alpha
		alpha1[i2] = L
		Lobj = DualobjFunc(alpha1,data)
		alpha1[i2]=H
		Hobj = DualobjFunc(alpha1,data)
		if Lobj > Hobj+eps:
			a2 = L
		elif Lobj < Hobj-eps:
			a2 = H
		else:
			a2 = alph2
	if abs(a2-alph2) < eps*(a2+alph2+eps):
		return 0
	a1 = alph1 + s*(alph2-a2)
	alpha[i1]=a1
	alpha[i2]=a2

	# updating threshold 
	b1 = b[0]-E1-y1*(a1-alph1)*k11-y2*(a2-alph2)*k12
	b2 = b[0]-E2-y1*(a1-alph1)*k12-y2*(a2-alph2)*k22
	if b1==b2:
		b[0]=b1
	else:
		if a1 < C and a1 > 0:
			b[0] = b1
		elif a2 < C and a2 > 0:
			b[0] = b2
		else:
			b[0] = (b1+b2)/2
	# end
	return 1

def secondHueristic(alpha,data,b,c,error,i2,output):
	m,n=data.shape
	error = findError(alpha,data,b,C,output)
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
	

def examineExample(i2,alpha,data,b,C,error,output,primal,dual):
	tol = 0.0001
	m,n = data.shape
	error = findError(alpha,data,b,C,output)
	alph2 = alpha[i2]
	y2 = data[i2][n-1]
	E2 = error[i2]
	r2 = E2*y2
	x = DualobjFunc(alpha,data)
	y = PrimalobjFunc(alpha,data,b,C,output)
	dual.append(x)
	primal.append(y)
	print x, y
	if (r2 < -tol and alph2 < C) or (r2 > tol and alph2 > 0):
		# finding no of non-zero, non-C alphas
		ind = np.where((alpha > 0) & (alpha < C))[0]
		indlist = ind.tolist()
		k = len(indlist)
		if k!=0:
			i1 = secondHueristic(alpha,data,b,C,error,i2,output)
			if takeStep(i1,i2,alpha,data,b,C,error,output):
				return 1
		
		if k!=0:
			for i1 in indlist:
				if takeStep(i1,i2,alpha,data,b,C,error,output):
					return 1
		
		for i1 in range(m):
			if i1!=i2:
				if takeStep(i1,i2,alpha,data,b,C,error,output):
					return 1
	return 0


def plotdata(data):
	m,n = data.shape
	data1_x=[]
	data1_y=[]
	data2_x=[]
	data2_y=[]
	plt.figure(1)
	for i in range(m):
		if data[i][n-1]==1:
			data1_x.append(data[i][0])
			data1_y.append(data[i][1])
		elif data[i][n-1]==-1:
			data2_x.append(data[i][0])
			data2_y.append(data[i][1])
	plt.plot(data1_x,data1_y,'ro')
	plt.plot(data2_x,data2_y,'bs')
	plt.xlim(-1.5,1.5)
	plt.ylim(-1.5,1.5)
	plt.xlabel('------> x')
	plt.ylabel('------> y')
	plt.show()

if __name__ == '__main__':
	data = np.array([[1,1,-1],[1,-1,1],[-1,-1,-1],[-1,1,1]])
	(m,n)=data.shape
	plotdata(data)
	alpha = np.zeros(m)
	error = np.zeros(m)
	output = np.zeros(m)
	b = np.zeros(1)
	C = 0.5
	max_iter = 50
	num_iter=0
	numChanged = 0
	examineAll = 1
	primal = []
	dual = []
	while numChanged > 0 or examineAll:
		numChanged = 0

		if examineAll:
			for i in range(m):
				numChanged += examineExample(i,alpha,data,b,C,error,output,primal,dual)
		else:
			for i in range(m):
				if alpha[i] > 0 and alpha[i] < C:
					numChanged += examineExample(i,alpha,data,b,C,error,output,primal,dual)
		if examineAll == 1:
			examineAll = 0
		elif numChanged == 0:
			examineAll = 1
		print output
		ind = np.where(output==1)[0]
		if num_iter > max_iter:
			break
	for i in range(m):
		y2 = data[i][n-1]
		point = data[i,0:n-1]
		output[i] = np.sign(SVMoutput(data,alpha,point,b)*y2)
	
	print "threshold: ", b
	print "Alpha: ", alpha
	print "output: ", output
	plt.figure(2)
	for i in range(m):
		if alpha[i]==0:
			if data[i][n-1]==1:
				plt.plot(data[i][0],data[i][1],'ro')
			else:
				plt.plot(data[i][0],data[i][1],'r+')
		else:
			if data[i][n-1]==1:
				plt.plot(data[i][0],data[i][1],'bo')
			else:
				plt.plot(data[i][0],data[i][1],'b+')
	
	# plotting decision boundary
	xmin = min(data[:,0])
	xmax = max(data[:,0])
	ymin = min(data[:,1])
	ymax = max(data[:,1])
	#print xmin, xmax, ymin, ymax
	boundary_x_0 = []
	boundary_y_0 = []
	boundary_x_1 = []
	boundary_y_1 = []
	boundary_x_2 = []
	boundary_y_2 = []
	for i in np.linspace(xmin-0.5, xmax+0.5, 400):
		for j in np.linspace (ymin-0.5, ymax+0.5, 400):
			point = np.array([i,j])
			out = SVMoutput(data,alpha,point,b)
			if out > -0.001 and out < 0.001:
				boundary_x_0.append(i)
				boundary_y_0.append(j)
			elif out > -1-0.001 and out < -1+0.001:
				boundary_x_1.append(i)
				boundary_y_1.append(j)
			elif out > 1-0.001 and out < 1+0.001:
				boundary_x_2.append(i)
				boundary_y_2.append(j)
	
	#k = min(len(boundary_x_0),len(boundary_x_1),len(boundary_x_2))
	#print k
	plt.plot(boundary_x_0[:],boundary_y_0[:],'g--') 
	#plt.plot(boundary_x_1[:k],boundary_y_1[:k],'b--')
	#plt.plot(boundary_x_2[:k],boundary_y_2[:k],'r--')
	plt.xlabel('-----------> x')
	plt.ylabel('-----------> y')
	plt.xlim(-1.5, 1.5)
	plt.ylim(-1.5, 1.5)
	plt.show()
	
	#plot dual objective value vs iterations
	'''
	k = min(len(primal), len(dual))
	plt.figure(3)
	plt.plot(range(k),primal[:k],'r-',label='Primal obj func')
	plt.plot(range(k),dual[:k],'b-',label='Dual obj func')
	plt.xlabel('------> iterations')
	plt.ylabel('------> objective function')
	plt.title('convergence plot')
	plt.axis([0,k,min(min(primal),min(dual)),max(max(primal),max(dual))])
	plt.show()
	'''

