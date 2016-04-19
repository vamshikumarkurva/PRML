import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm

if __name__=='__main__':
	data = np.loadtxt('/home/iistlab/Desktop/vamshi/CASP1.csv',delimiter=",")
	m,n=data.shape
	X=data[:,0:n-1]
	Y=data[:,n-1]
	clf = svm.SVR(gamma=150,kernel='rbf')
	clf.fit(X,Y)
	E=np.zeros(m)
	if n==2:
		x=[]
		y=[]
		z=[]
		for i in range(m):
			x.append(X[i].item())
			k = clf.predict(X[i]).item()
			err = abs(k-Y[i])
			E[i] = err
			#print err
			y.append(k)
			z.append(Y[i].item())
		print "Mean square error: ", E.std()
		plt.plot(x,y,'rs',label='predicted')
		plt.plot(x,z,'bo',label='original')
		plt.show()
	else:				
		for i in range(m):
			k = clf.predict(X[i]).item()
			err = abs(k-Y[i])
			E[i] = err
			print err
		print "Mean square error: ", E.std()
