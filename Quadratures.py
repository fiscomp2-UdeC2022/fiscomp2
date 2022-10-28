import numpy as np

__all__ = ["Rectangle", "MidPoint", "Trapezoide", "Romberg"]

def Rectangle(f, a, b, N=50):
    x = np.linspace(a,b,N+1)
    h = x[1] - x[0]
    
    return h * np.sum( f(x[:-1]) ) 


def MidPoint(f, a, b, N=50):
    x = np.linspace(a,b,N+1)
    h = x[1] - x[0]
    
    return h * np.sum( f(x[:-1] + h/2) ) 


def Trapezoide(f, a, b, N=50):
    x = np.linspace(a,b,N+1)
    h = x[1] - x[0]

    fval = f(x)

    return h*fval.sum() - 0.5*h*(fval[0]+fval[-1])


def _Richardson(Q1, Q0, m):
    """Extrapola Q1 y Q0 usando extrapolación de Richardson de orden m"""
    k = 4.0**m

    return (k * Q1 - Q0)/(k-1.0)


def Romberg(f, a, b, tol=1e-8, rtol=1e-8, maxiter=10, show=False):
    h = b-a
    
    N = 1
    Q0 = [Trapezoide(f, a, b, N)]

    if show:
        print(f"iter 0: {h}\t {Q0[0]}")
        
    for i in range(1, maxiter+1):
        N = N*2
        h = 0.5*h

        x = a + np.arange(1,N,2) * h
        Q1 = [0.5 * Q0[0] + h * np.sum( f(x) )]

        for j in range(i):
            Q1 += [ _Richardson(Q1[j], Q0[j], j+1) ]

        if show:
            print(f"iter {i}: {h}\t {Q1}")

        result = Q1[-1]
        error = np.abs(result - Q0[-1])
        if error<tol or error < rtol*np.abs(result):
            break
        
        Q0 = Q1
        
    else:
        print(f"se excede el maximo numero de iteraciones {maxiter}")

    return result

