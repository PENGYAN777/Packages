#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:39:21 2024

bisection method for a single equation

@author: user
"""
import CoolProp as CP
import time
import os
from IPython import get_ipython;   
get_ipython().magic('reset -sf')
os.system('clear')
fluidname = "HEOS::MD4M"
start = time.time()

def bisection_method(f,p1,p2,p3, a, b, tol=1e-4, max_iterations=100):
    """
    Bisection method to find a root of the function f(x) in the interval [a, b].

    Parameters:
        f: Callable function representing the equation to find the root of.
        a: Lower bound of the interval.
        b: Upper bound of the interval.
        tol: Tolerance level for the root (default is 1e-6).
        max_iterations: Maximum number of iterations allowed (default is 1000).

    Returns:
        root: Approximate root of the function within the specified tolerance.
        iterations: Number of iterations taken to find the root.
    """
    if f(a,p1,p2,p3) * f(b,p1,p2,p3) >= 0:
        # print("a,b ", a,b)
        # print("p1,p2,p3 ", p1,p2,p3)
        print("f(a), f(b) ", f(a,p1,p2,p3),f(b,p1,p2,p3) )
        raise ValueError("The function values at the endpoints must have opposite signs.")

    iteration = 0
    while iteration < max_iterations:
        c = (a + b) / 2  # Calculate midpoint
        if abs(f(c,p1,p2,p3)) < tol:
            return c, iteration  # Root found within tolerance
        elif f(a,p1,p2,p3) * f(c,p1,p2,p3) < 0:
            b = c  # Update upper bound
        else:
            a = c  # Update lower bound
        iteration += 1

    raise ValueError("Bisection method did not converge within the maximum number of iterations.")


# Example usage:
# if __name__ == "__main__":
#     import math

#     # Define the function for which we want to find the root
#     def func(x,p1,p2,p3):
#         return x**3 - 2*x - p1-p2-p3

#     # Initial interval
#     a = 1
#     b = 3
#     p1 = 6
#     p2 = 6
#     p3 = 6

#     # Find the root using the bisection method
#     root, iterations = bisection_method(func, p1,p2,p3, a, b)

#     # Print the result
#     print("Approximate root:", root)
#     print("Number of iterations:", iterations)

end = time.time()
print("computational time(s): ", end - start)
