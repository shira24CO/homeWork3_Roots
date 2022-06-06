import sys

import mpmath
import mpmath.math2
import sympy as sp
from sympy.utilities.lambdify import lambdify

epsilon = 0.000001
#print the derivative in a symbol way
x =sp.symbols('x')
my_f =x**4 +x**3 -3*x**2
print("my_func: ", my_f)
my_f1=sp.diff(my_f, x)
print("f' : ", my_f1)

d1 = sp.diff(my_f1,x)
print("f'': ", d1)

#calc the derivative from func -> lambdify
f = x**4 +x**3 -3*x**2
f_prime = f.diff(x)
print("f : ",f)
print("f' : ",f_prime)
f = lambdify(x, f)
f_prime = lambdify(x, f_prime)
d1 = lambdify(x,d1)
print("f(1):",f(1))
print("f'(1):",f_prime(1))
print("f''(1)",d1(1))

def bisection(f, start_point, end_point):
    a_n = start_point
    b_n = end_point
    counter = 0
    while abs(b_n - a_n) >= epsilon and counter < int((-1)*(mpmath.ln(epsilon/(end-start)))/mpmath.ln(2)-1):
        m_n = (a_n + b_n) / 2
        f_m_n = f(m_n)
        counter += 1
        if f(a_n) * f_m_n < 0:
            # a_n = a_n  # אפשר למחוק?
            b_n = m_n
        elif f(b_n) * f_m_n < 0:
            a_n = m_n
            # b_n = b_n # אפשר למחוק?
        elif f_m_n == 0:
            print("Found exact solution.")
            return m_n
        else:
            print("Bisection method fails.")
            return None
    print('Found solution after', counter, 'iterations.')
    # return (a_n + b_n) / 2
    return a_n



def newton(f, Df, x0, max_iter):
    xn = x0
    iteration_counter = 0
    while abs(xn-f(xn)) > epsilon and iteration_counter < 100:
        fxn = f(xn)
        if abs(fxn) < epsilon:
            print('Found solution after', iteration_counter, 'iterations.')
            return xn
        Dfxn = Df(xn)
        if Dfxn == 0:
            print('Zero derivative. No solution found.')
            return None
        iteration_counter += 1
        xn = xn - (fxn / Dfxn)
    print('Exceeded maximum iterations. No solution found.')
    return None


def secant(f, x0, x1):
    f_x0 = f(x0)
    f_x1 = f(x1)
    iteration_counter = 0
    while abs(f_x1) > epsilon and iteration_counter < 100:
        try:
            denominator = float(f_x1 - f_x0) / (x1 - x0)
            x = x1 - float(f_x1) / denominator
        except ZeroDivisionError:
            print( "Error! - denominator zero for x = ", x)
            sys.exit(1)  # Abort with error
        x0 = x1
        x1 = x
        f_x0 = f_x1
        f_x1 = f(x1)
        iteration_counter += 1
    # Here, either a solution is found, or too many iterations
    if abs(f_x1) > epsilon:
        iteration_counter = -1
    if iteration_counter > 0:  # Solution found
        print("Number of function calls: %d" % (2 + iteration_counter))
        return x
    else:
        print("Solution not found!")



# switch_case
def switch_case(choise, start, end):
    i = start
    #1 for Bisection,2 for  Newton Rapson , 3 for secant
    if choise == '1':
        if f(0) == 0:
            print("A solution: ",f(0))
        while i < end:
            if f(i) * f(i + 0.1) < 0:
                print(bisection(f, i, i+0.1))
            i += 0.1
        i = start
        value = 0
        while i < end:
            if f_prime(i) * f_prime(i+0.1) < 0:
                value = bisection(f_prime,i,i+0.1)
                if abs(f(value)) < epsilon:
                    print(value,"is ~ zero point", f(value))
                else:
                    print("The derivative is zero, but the value of the function is different than zero")
                #print(bisection(f_prime, i, i+0.1))
                    #print(bisection(f_prime,i,i+0.1))
            i = i + 0.1

    elif choise == '2':
        if f(0) == 0:
            print("A solution: ",f(0))
        while i < end:
            if f(i) * f(i + 0.1) < 0:
                print(newton(f, f_prime, i, 100))
            i = i + 0.10
        i = start
        value = 0
        while i < end:
            if f_prime(i) * f_prime(i + 0.1) < 0:
                value = newton(f_prime,d1,i,100)
                if value is None:
                    exit(1)
                if abs(f_prime(value)) < epsilon:
                    print(value,"is ~ zero point", f_prime(value))
                else:
                    print("The derivative is zero, but the value of the function is different than zero")
                # print(bisection(f_prime, i, i+0.1))
                # print(bisection(f_prime,i,i+0.1))
            i = i + 0.1
    elif choise == '3':
        if f(0) == 0:
            print("A solution: ",f(0))
        while i < end:
            if f(i) * f(i + 0.1) < 0:
                solution = secant(f, i, i + 0.1)
                print("A solution is: %f" % (solution))
            i = i + 0.1
        i = start
        while i < end:
            if f_prime(i) * f_prime(i + 0.1) < 0:
                solution = secant(f_prime, i, i + 0.1)
                if solution is None:
                    exit(1)
                if abs(f(solution)) < 1.0e-6:
                    print("A solution prime is: %f" % (solution))
                    print(solution, "is ~ zero point", f(solution))
            i = i + 0.1

    else:
        print("Bad Choise")


start = input("enter  start range: ")
start = int(start)
end = input("enter end range: ")
end = int(end)
choise = input(
    "In what method you interested in finding the roots of the equation \n1 for Bisection, 2 for  Newton Rapson , 3 for secant\nyour choise ? : ")
switch_case(choise,start,end)