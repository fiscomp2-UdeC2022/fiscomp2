import numpy as np
from GaussianElimination import Resolver


class simple_cspline:
    def __init__(self, xi, yi):
        """Calcula los coeficientes de polinomios cubicos `pi(x) = ai + bi
        (x-xi) + ci (x-xi)^2 + di (x-xi)^3` entre dos valores
        consecutivos de los nodos `xi` asumiendo que `pi(xi) =
        yi`. Esta rutina estima la derivada `pi'(xi) = dyi` usando
        derivadas centradas y suponiendo que la derivada es cero en los extremos `dy0=0=dyN`.

        """
        
        self.xi = np.asarray(xi)
        yi = np.asarray(yi)

        
        dy = (yi[2:] - yi[:-2])/(self.xi[2:] - self.xi[:-2])
        dy = np.concatenate([[0], dy, [0]])

        dx = self.xi[1:] - self.xi[:-1]
        A = ( (yi[1:]-yi[:-1])/dx - dy[:-1]) / dx
        B = (dy[1:] - dy[:-1])/dx
        
        self.coef = np.array([yi[:-1],
                              dy[:-1],
                              3*A - B,
                              (B-2*A)/dx]).T
        
    
    def __call__(self, x):
        x = np.asarray(x)

        i = np.where(
            (self.xi[:-1,None] <= x[None,:])*(self.xi[1:,None]  >  x[None,:]))[0]

        t = x - self.xi[i]

        a, b, c, d = self.coef[i].T
        return a + t*(b + t*(c + d*t))


class cspline:
    """Calcula los coeficientes de polinomios cubicos `pi(x) = ai + bi
    (x-xi) + ci (x-xi)^2 + di (x-xi)^3` entre dos valores consecutivos
    de los nodos `xi` asumiendo que `pi(xi) = yi`. Esta rutina calcula
    los coeficientes asumiendo que tanto `pi(x)` como `pi'(x)` son
    continuos en los nodos.

    """
    def __init__(self, xi, yi):
        self.xi = np.asarray(xi)
        yi = np.asarray(yi)
        
        pendiente = np.empty(yi.shape)

        pendiente[0]    = 3*(yi[1] - yi[0])
        pendiente[1:-1] = 3*(yi[2:] - yi[:-2] )                                     
        pendiente[-1]   = 3*(yi[-1] - yi[-2])

        matriz = np.zeros((self.xi.size, self.xi.size))

        matriz[0,:2] = [2.0, 1.0]
        matriz[-1, -2:] = [1.0, 2.0]

        for i in range(1, self.xi.size-1):
            matriz[i, i-1:i+2] = [1.0, 4.0, 1.0]

        a = Resolver(matriz, pendiente)

        b = 3*(yi[1:] - yi[:-1]) - 2*a[:-1] - a[1:]
        c = a[:-1] + a[1:] - 2*(yi[1:] - yi[:-1])

        self.coef = np.array([yi[:-1], a[:-1], b, c]).T
        

    def __call__(self, x):
        x = np.asarray(x)

        i = np.where(
            (self.xi[:-1,None] <= x[None,:])*(self.xi[1:,None]  >  x[None,:]))[0]

        t = (x - self.xi[i]) / (self.xi[i+1] - self.xi[i]) 

        a, b, c, d = self.coef[i].T

        return a + t*(b + t*(c + d*t))


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # nodos para la interpolacion
    xi = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    yi = [1.13e-7, 1.24e-4, 0.0183, 0.368, 1.0, 0.368, 0.0183, 1.23e-4, 1.13e-7]


    # puntos donde evaluaremos la funcion interpolada
    x = np.linspace(np.min(xi), np.max(xi), 200, endpoint=False)
    
    # calcula los coeficientes de los polinomios entre dos nodos usando dos metodos distintos
    simple_sp = simple_cspline(xi, yi) 
    csp = cspline(xi, yi)

    print("coeficientes usando `cspline`:\n", csp.coef, "\n")
    print("coeficientes usando `simple_cspline`:\n", simple_sp.coef)

    # las variables `simple_sp` y `csp` son funciones que pueden ser evaluadas en cada `x`
    plt.plot(x, simple_sp(x), label='simple_cspline')
    plt.plot(x, csp(x), label='cspline tradicional')

    # grafica los nodos
    plt.scatter(xi, yi, label="nodos")


    plt.legend()
    plt.show()

    
