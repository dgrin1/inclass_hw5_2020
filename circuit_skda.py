from numpy import array,empty,copy,dot
from numpy.linalg import inv,solve
#scipy.linalg lu
#dot

A=array([[1,1,-1],
	 [-6,0,0],
	 [0,0,-9]],float)
v=array([0,6,-12],float)
N=len(v)
b=copy(A)
c=copy(v)

# Gaussian elimination
for m in range(N):

    # Divide by the diagonal element
    div = A[m,m]
    A[m,:] /= div
    v[m] /= div

    # Now subtract from the lower rows
    for i in range(m+1,N):
        mult = A[i,m]
        A[i,:] -= mult*A[m,:]
        v[i] -= mult*v[m]

# Backsubstitution
x = empty(N,float)
for m in range(N-1,-1,-1):
    x[m] = v[m]
    for i in range(m+1,N):
        x[m] -= A[m,i]*x[i]

print(x)



import scipy.linalg as la
# Get LU decomposition and save the L & U matrices, where p is ???
p,l,u=la.lu(b)

#Compute solution to Ly=c
y=dot(inv(l),dot(p,c))
#Compute solution to Ux=y
newx=dot(inv(u),y)
print(newx)








