__all__ = ["Vandermonde", "Lagrange", "Neville"]

import numpy as np
from GaussianElimination import Resolver
from PolyClass import Poly


class Vandermonde(Poly):
    def __init__(self, xi, yi):
        """Encuentra los coeficientes a[n] de un polinomio p(x)=sum( a[n]*x^n,
        n=0,N), con N el número de puntos {xi,yi}."""

        xi = np.asarray(xi)
        yi = np.asarray(yi)
        
        self.dim = xi.size
        
        # crea una matriz de Vandermonde
        MatrizVandermonde = np.vander(xi, increasing=True)

        # usa eliminacion gaussiana para encontrar los coeficientes a[n]
        # luego, guarda cada coeficiente en una lista
        self.coef = Resolver(MatrizVandermonde, yi)

        
class Lagrange(Poly):
    def __init__(self, xp, yp):
        """Clase para interpolar usando método de Lagrange"""
        super().__init__([0])
        
        self.xp = np.asarray(xp)
        self.yp = np.asarray(yp)

        L = [Poly([-xj,1]) for xj in xp]

        pol = Poly([0.0])

        # para cada valor de i, calcula polinomios de Lagrange
        for i,(xi,yi) in enumerate(zip(xp,yp)):
            tmp = np.prod(L[:i]) * np.prod(L[i+1:])
            tmp = tmp * (yi / tmp(xi))

            pol = pol + tmp

        self.coef = pol.coef

        
class Neville(Poly):
    def __init__(self, xp, yp):
        super().__init__([0])
        
        self.xp = np.asarray(xp)
        self.yp = np.asarray(yp)

        p = [Poly([yi]) for yi in yp]
        
        for j in range(1,len(yp)):
            for i in range(len(yp)-j):
                factor = 1.0/(xp[i]-xp[i+j])
                pol1 = Poly([-xp[i+j], 1])
                pol2 = Poly([-xp[i], 1])

                p[i] = factor*( pol1*p[i] - pol2*p[i+1])

        self.coef = p[0].coef




# esto es para testear el código directamente
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
    # xp = np.array([1, 2, 5, 7, 8, 9], float)
    # yp = np.array([1, -2, -1, -2, 0, 2], float)
    npoints = 20
    xp = np.linspace(1, 9, npoints)
    yp = np.random.randn(npoints)


    for metodo in [Lagrange, Neville, Vandermonde]:
        p = metodo(xp, yp)
    
        print(f"Polinomio interpolado: \ngrado={p.deg}\ncoef={p}")

        x = np.linspace(0.5, 9.5, 500)
        plt.plot(x, p(x), label=f"metodo de {metodo.__name__}")

    plt.scatter(xp, yp, label="puntos interpolados")
    plt.legend()
    plt.title("Interpolación de Lagrange")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.ylim(-3,3)
    plt.show()


        
