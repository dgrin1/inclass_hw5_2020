from numpy import exp

guess = 1

def f(x):
    return float(exp(1-x**2))

while (abs(guess-f(guess)) > 0.01):
    guess = f(guess)

print(guess)