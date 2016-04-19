import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm

if __name__=='__main__':
	data = np.loadtxt('/home/iistlab/Desktop/vamshi/assign2/data2.csv',delimiter=",")
	m,n=data.shape
	X=data[:,0:n-1]
	Y=data[:,n-1]
	clf = svm.SVC(gamma=100,kernel='rbf')
	clf.fit(X,Y)
	h = 0.02 # step size
	xmin=min(X[:,0])-1
	xmax=max(X[:,0])+1
	ymin=min(X[:,1])-1
	ymax=max(X[:,1])+1
	xx,yy = np.meshgrid(np.arange(xmin,xmax,h),np.arange(ymin,ymax,h))
	Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
	Z = Z.reshape(xx.shape)
	plt.contourf(xx,yy,Z,alpha=0.2)
	for i in range(m):
		if Y[i]==1:
			plt.plot(X[i][0],X[i][1],'r*')
		elif Y[i]==-1:
			plt.plot(X[i][0],X[i][1],'bo')
	plt.xlabel('---------> X')
	plt.ylabel('---------> X')
	plt.xlim(xx.min(),xx.max())
	plt.ylim(yy.min(),yy.max())
	plt.xticks(())
	plt.yticks(())
	plt.title('Decision Boundary')
	plt.show()
