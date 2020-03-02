# 02/20/2019   
# inclass work for integration -- Romberg Integral
import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate

# Trapezoid Sum (using equation in chapter 5 p142 of Newman)
def TrapezoidSum(f,a,b,N):
    delta=float(abs(b-a)/N)
    sum=0.0
    for i in range(1,N):
        # taking care of the case if a>b
        if b > a:
            sum+=f(a+i*delta)
        else:
            sum+=f(b+i*delta)
    sum=sum+(1/2)*(f(a)+f(b))
    sum*=delta
    # taking care of the case if a>b
    if b > a:
        return sum
    else:
        return sum*(-1.0)

# Romberg integation implemented according to chapter 5.4 p160 of Newman
def RombergIntegral(f,a,b,N):
    # Define a dictionary for storing the estimated integral R_(i,m)
    R={} 
    # Define the first two estimates of the integral using trapezoidal rule
    R[1,1]=TrapezoidSum(f,a,b,N)
    R[2,1]=TrapezoidSum(f,a,b,2*N)
    # Define a variable to update the iteration number for TrapezoidSum()
    Nmin=2*N
    # Define variables for keeping track of the indexes of the estimated Ri,m
    i=2
    m=1
    epsilon=1e-16 # target accuracy
    error=1e0 # variable for checking the computational error
    while abs(error) > epsilon: 
        # Conditional case that calculates R_i,m whenever i is larger than m
        if m<i:    
            # Defind the variable R_i,m+1 for iteration
            R[i,m+1]=R[i,m]+((1.0)/((4**m)-1))*(R[i,m]-R[i-1,m])
            # Compute the error using eqn-5.49 from chapter 5.4 p161 of Newman
            error = ((1.0)/((4**m)-1))*(R[i,m]-R[i-1,m])
            m+=1
        # Conditional case that calculates R_i,m whenever m is equal to i and i is upgraded 1 higher
        else:
            # Upgrade the value of the parameter i and reset the value of the parameter m
            i+=1
            m=1
            # Update the minimum iteration number for TrapezoidSum()
            Nmin*=2
            # Update the elements in the dictionary R_(i,m)
            R[i,m]=TrapezoidSum(f,a,b,Nmin)
    return R[i,m],error,i,Nmin/N

# test
'''
print("Uisng Trapezoid Integral with 10000 iterations gives the result:", TrapezoidSum(lambda x: x**4-2*x+1,0,2,10000))
print("Using Romberg Integral gives the result:", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[0],"with a trapezoidal integration level of", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[2],"and an iteration number of",RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[3])
print("The error using Romberg Integral is:", RombergIntegral(lambda x: x**4-2*x+1,0,2,10)[1])
'''

# Using the Romberg Integral to do thermal integral practice found in "thermal.py" or "thermal_gauss.py" -- details about the question see lecture note 02/19/2020 p7

# Define a function that performs integration over the infinite ranges of energy (Eb<=E<=inf) -- A=int(E_B,inf,(Eb)^2/(exp(Eb/T)-1)) with respect to Eb=13.6 eV
# Using change of variable trick shows in Chapter 5.8 p179 of Newman to define the integrand after changing the variable as a function of z (the changed variable of E) & T.
# Define the integrands for calculating the numerator A as well as the normalizing factor N such that Probability=A/N

# Integrand for A using the change of variable of x=z/(1-z)+a
def A_integrand1(T):
    # Define the constant Egamma as Eb as the lower limit for the original integral
    Eb=13.6 # eV
    # Returning the integrand as functions of z & T for both the integral A and the normalizing factor N
    return lambda z: (1/((1-z)**2))*((z/(1-z)+Eb)**2)/(np.exp((z/(1-z)+Eb)/T)-1)

# Integrand for N using the change of variable of x=z/(1-z)
def N_integrand1(T):
    return lambda z: (1/((1-z)**2))*((z/(1-z))**2)/(np.exp((z/(1-z))/T)-1)

# Using Romberg Integration and the integrand defined above to compute the probability
def Probability1(T):
    # Define the integration limits for the transformed integral -- note that b could not be 1 because of the divide by 0 error and a could not be 0 to avoid nan
    a=0.0000000001
    b=0.9999999999
    # Define the number of iterations for performing TrapezoidSum
    N=100
    # Calculate the integral and the error with the given T as the numerator for getting the probability (from Eb to inf) using Romberg Integral
    AResult = RombergIntegral(A_integrand1(T),a,b,N)
    A = AResult[0]
    Aerr = AResult[1]
    # Calculate the integral and the error with the given T as the denominator for getting the probability (Normalizing factor) using Romberg Integral
    NResult = RombergIntegral(N_integrand1(T),a,b,N)
    N = NResult[0]
    Nerr = NResult[1]
    # Caculate the probability using P=A/N
    P=A/N
    # Error calculated using error propagation -- might be a bad idea because the error is too small
    Perr1 = np.abs(P)*np.sqrt((Aerr/A)**2+(Nerr/N)**2)
    # Error calculated by simply average the two errors Aerr and Nerr
    Perr2 = (np.abs(Aerr)+np.abs(Nerr))/2
    # Error calculated based on the deviation from the original values gives only 0.0
    #Perr=float(np.abs(P-(A+Aerr)/(N-Nerr)))
    return P, Perr1, Perr2


# Test the result -- obtain the same result of probability when T=8.0 K as thermal_gausspy 
print("")
print("Thus, the probability when T=8.0 K is:", Probability1(8.0)[0], "with an estimated error of (1)",Probability1(8.0)[1],"using error propagation for division OR (2)",Probability1(8.0)[2],"using average error.")
print("")
print("Thus, the probability when T=10.0 K is:", Probability1(10.0)[0], "with estimated errors of",Probability1(10.0)[1],"using error propagation for division OR (2)",Probability1(10.0)[2],"using average error.")
print("")
print("Thus, the probability when T=50.0 K is:", Probability1(50.0)[0], "with estimated errors of",Probability1(50.0)[1],"using error propagation for division OR (2)",Probability1(50.0)[2],"using average error.")
print("")
# i.e. Thus, the probability when T=100.0 K is: 0.9963246818967767 with estimated errors of (1) 2.8037593259440285e-23 using error propogation for division between A & N; (2) 3.3827107907563535e-17 using the average of the errors from A & N.
print("For some reasons, the higher the value of T is, the longer it takes for Romberg Integration to compute the result. ")
print("")


# Define a function that plots the probability for a range of temperature
def PlotP1():
    # Initialize the variable T from 1 to 1000.5 (same as the karr of thermal.py) as an array
    T=np.arange(1.0,10.6,0.05)
    # Initialize the variable Result for storing the results for probability and errors
    Result=np.array([])
    for Temp in T:
        # Calcualting the probability & error -- for each iteration, only runs the function Probability(T) once for optimization
        Result = np.append(Result,Probability1(Temp))
    # Selecting and grouping the data points and storing them as different variables (i.e. Probability and Error)
    size = np.size(Result)
    P = Result[0:size:3]
    Perr = Result[2:size:3]
    
    # Clear the entire current figure in case of multiple runnings of the function
    plt.clf() 
    # Enable interactive mode
    plt.ion()
    # Use Latex font serif
    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    
    ## Two panel subplots
    fig, ax = plt.subplots(2,sharex=True)
    # Plot "Probability vs. Temperature" on subplot 1
    ax[0].plot(T,P,'-',color='#ffaf0f80',label='Temperature-Dependent \n Probability of ionization',linewidth=2)
    # Subplot axis labels
    ax[0].set_ylabel('Probability of Ionization in Log scale')
    ax[0].set_yscale('log')

    # Plot the average error on subplot 2
    ax[1].plot(T,Perr,'-',color='k',label='Probability Error',linewidth=2)
    # Subplot axis labels
    ax[1].set_ylabel('Errors for Probability')

    # Enable grid
    ax[0].grid(linestyle='-', linewidth=0.1)
    ax[1].grid(linestyle='-', linewidth=0.1)
    # Include legends separately for suplots
    ax[0].legend(loc=4)
    ax[1].legend(loc=1)
    
    # Common x label
    plt.xlabel('Temperatue (K)')
    # Title for both subplots
    fig.suptitle('Probability of ionization with errors vs. Temperature')
    
    # Show & save the figure
    fig.show()
    fig.savefig('HYL_thermal_Romberg_integral1.pdf',format='pdf')


# Plot
PlotP1()





'''
Try another kind of integrand (trignometric) to see if the code will be faster
'''

'''

# Integrand for A using the trignometric change of variable of 
def A_integrand2(T):
    
    # Define the constant Egamma as Eb as the lower limit for the original integral
    Eb=13.6 # eV
    
    # Returning the integrand as functions of z & T for both the integral A and the normalizing factor N
    # return lambda theta: (1/(np.cos(theta)**2))*((np.tan(theta)+Eb)**2)/(np.exp(np.tan(theta)+Eb)/T-1)
    
    return lambda theta: ((Eb+np.tan(theta))**2.0)*(np.exp(-np.tan(theta)))/(np.exp(Eb)-np.exp(-np.tan(theta)))*(np.power(np.cos(theta),-2.e0))

# Integrand for N using the trignometric change of variable of 
def N_integrand2(T):
    return lambda theta: ((0.0+np.tan(theta))**2.0)*(np.exp(-np.tan(theta)))/(np.exp(0.0)-np.exp(-np.tan(theta)))*(np.power(np.cos(theta),-2.e0))

# Using Romberg Integration and the trignometirc integrand defined above to compute the probability
def Probability2(T):
    # Define the integration limits for the transformed integral -- note that b could not be 1 because of the divide by 0 error and a could not be 0 to avoid nan
    a=0.000000001
    b=np.pi/2-a
    # Define the number of iterations for performing TrapezoidSum
    N=100
    # Calculate the integral and the error with the given T as the numerator for getting the probability (from Eb to inf) using Romberg Integral
    AResult = RombergIntegral(A_integrand2(T),a,b,N)
    A = AResult[0]
    Aerr = AResult[1]
    # Calculate the integral and the error with the given T as the denominator for getting the probability (Normalizing factor) using Romberg Integral
    NResult = RombergIntegral(N_integrand2(T),a,b,N)
    N = NResult[0]
    Nerr = NResult[1]
    # Caculate the probability using P=A/N
    P=A/N
    # Error calculated using error propagation -- might be a bad idea because the error is too small
    Perr1 = np.abs(P)*np.sqrt((Aerr/A)**2+(Nerr/N)**2)
    # Error calculated by simply average the two errors Aerr and Nerr
    Perr2 = (np.abs(Aerr)+np.abs(Nerr))/2

    # Error calculated based on the deviation from the original values gives only 0
    #Perr=float(np.abs(P-(A+Aerr)/(N-Nerr)))
    return P, Perr1, Perr2


# Test the result -- obtain the same result of probability when T=8.0 K as thermal_gausspy 
print("")
print("Thus, the probability when T=8.0 K is:", Probability2(100.0)[0], "with an estimated error of",Probability2(100.0)[1],"using error propagation for division OR (2)",Probability2(8.0)[2],"using average error.")


# Define a function that plots the probability for a range of temperature
def PlotP2():
    # Initialize the variable T from 1 to 1000.5 (same as the karr of thermal.py) as an array
    T=np.arange(1.0,10.6,0.05)
    # Initialize the variable Result for storing the results for probability and errors
    Result=np.array([])
    for Temp in T:
        # Calcualting the probability & error -- for each iteration, only runs the function Probability(T) once for optimization
        Result = np.append(Result,Probability2(Temp))
    # Selecting and grouping the data points and storing them as different variables (i.e. Probability and Error)
    size = np.size(Result)
    P = Result[0:size:3]
    Perr = Result[2:size:3]
    
    # Clear the entire current figure in case of multiple runnings of the function
    plt.clf() 
    # Enable interactive mode
    plt.ion()
    # Use Latex font serif
    plt.rc('text',usetex=True)
    plt.rc('font', family='serif')
    
    ## Two panel subplots
    fig, ax = plt.subplots(2,sharex=True)
    # Plot "Probability vs. Temperature" on subplot 1
    ax[0].plot(T,P,'-',color='#ffaf0f80',label='Temperature-Dependent \n Probability of ionization',linewidth=2)
    # Subplot axis labels
    ax[0].set_ylabel('Probability of Ionization in Log scale')
    ax[0].set_yscale('log')

    # Plot the average error on subplot 2
    ax[1].plot(T,Perr,'-',color='k',label='Probability Error',linewidth=2)
    # Subplot axis labels
    ax[1].set_ylabel('Errors for Probability')

    # Enable grid
    ax[0].grid(linestyle='-', linewidth=0.1)
    ax[1].grid(linestyle='-', linewidth=0.1)
    # Include legends separately for suplots
    ax[0].legend(loc=4)
    ax[1].legend(loc=1)
    
    # Common x label
    plt.xlabel('Temperatue (K)')
    # Title for both subplots
    fig.suptitle('Probability of ionization with errors vs. Temperature')
    
    # Show & save the figure
    fig.show()
    fig.savefig('HYL_thermal_Romberg_integral2.pdf',format='pdf')

'''