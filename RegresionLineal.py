from numpy import vander, sum
import PolyClass
from GaussianElimination import Resolver

class polyfit(PolyClass.Poly):
    def __init__(self, x, y, grado=1):
        "Esta rutina realiza regresion lineal entre los puntos (x,y) y un polinomio p(x)=a[0]+a[1]*x+...+a[N]*x^N con N=grado, tal que la norma chi=sum( |yi-p(xi)|^2)  es minima. El metodo guarda los coeficientes a[n]."

        # Crea una matriz de Vandermonde. Esta matriz no es necesariamente cuadrada, sus dimensiones son (len(x),grado+1)
        G = vander(x, grado + 1, increasing=True)

        M = G.T @ G
        z = G.T @ y

        # encuentra los coeficientes del polinomio que minimiza el error con los datos (x,y)
        self.coef = Resolver(M, z)

        # evalua el error cuadratico
        self.chi2 = sum(( G @ self.coef - y )**2)
        
