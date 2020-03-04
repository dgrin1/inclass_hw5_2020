import numpy as np

# Include and upgrade the function of doing Gaussian elimination & back-substitution with partial pivoting for a linear system Ax=v
def gausselim(AA,vv):
    # Making copies so that the input arrays are not modified
    A=np.copy(AA)
    v=np.copy(vv)
    N = len(v)
    # Gaussian elimination
    for m in range(N):
        # Divide by the diagonal element
        div = A[m,m]
        # Check to see if partial pivoting is required
        if np.abs(div)<np.max(np.abs(A[m:,m])):
            SwapRow=np.where(np.abs(A[:,m])==np.max(np.abs(A[m:,m])))[0][0]
            A[[m,SwapRow]]=A[[SwapRow,m]]
            v[[m,SwapRow]]=v[[SwapRow,m]]
        # Reassign the diagonal element with a different value after partial pivoting
        div = A[m,m]
        A[m,:] /= div
        v[m] /= div
        # Now subtract from the lower rows
        for i in range(m+1,N):
            mult = A[i,m]
            A[i,:] -= mult*A[m,:]
            v[i] -= mult*v[m]
    # Backsubstitution
    x = np.empty(N,float)
    for m in range(N-1,-1,-1):
        x[m] = v[m]
        for i in range(m+1,N):
            x[m] -= A[m,i]*x[i]
    return x

'''
Test using Exercise 6.4 & inclass work
'''
A=np.array([[4,-1,-1,-1],
           [-1,3,0,-1],
           [-1,0,3,-1],
           [-1,-1,-1,4]],float)
# In the order of V1,V2,V3,V4
v = np.array([5, 0, 5, 0],float)
# Testing the updated version of gausselim(A,v) (including partial pivoting) using Exercise 6.4
print(np.linalg.solve(A,v))
print(np.dot(A,np.linalg.solve(A,v)))
print(gausselim(A,v))

print("")

B = np.array([[ 2,  1,  4,  1 ],
            [ 3,  4, -1, -1 ],
            [ 1, -4,  1,  5 ],
            [ 2, -2,  1,  3 ]], float)
b = np.array([ -4, 3, 9, 7 ],float)
# The correct answer should be [2,-1,-2,1] without undoing the swapping of the solution
print(np.linalg.solve(B,b))
print(np.dot(B,np.linalg.solve(B,b)))
# Solve the linear system with our gausselim()
print(gausselim(B,b))
