import numpy as np

class Poly:
    """Clase de polinomios. 
Para un polinomio de la forma `P(x)=p0 + p1*x + p2*x**2 + ... + pN*x**N`, esta clase guarda los coeficientes `coef=[p0 p1 ... pN]` y el grado del polinomio `deg=N`.

[Atributos]
    coef : Coeficientes del polinomio.
    deg  : Grado del polinomio.

[Metodos]
    __call__(x) : Retorna el valor del polinomio evaluado en `x`.

    deriv(m=1)  : Retorna la m-ésima derivada del polinomio.
"""
    def __init__(self, coef):
        self.coef = coef        # llama coef.setter

    @property
    def coef(self):
        """Retorna  """
        return self._coef      

    @coef.setter
    def coef(self, coef):
        self._coef = np.asarray(coef)
        self.deg = self.coef.size - 1

    def __call__(self, x):
        """La clase es tratada como un functor o llamable (callable). Si `p` es una instancia de `Poly`, entonces `p(x)` evalúa el polinomio en un valor de `x` usando el método de Horner. """
        x = np.asarray(x)

        value = np.zeros_like(x)   #self.coef[-1]

        for c in self.coef[::-1]:
            value = c + value * x

        return value

    
    def __pos__(self):
        return Poly(self.coef)

    
    def __neg__(self):
        coef = -self.coef
        return Poly(coef)

    
    def __add__(self, other):
        """Retorna la suma de dos polinomios. Si `other` es un número, se trata como un polinomio de grado cero `Poly([other])`"""
        if isinstance(other, (int, float, complex)):
            other = Poly([other])

        if self.deg < other.deg:
            pmin, pmax = self.coef, other.coef.copy()
        else:
            pmin, pmax = other.coef, self.coef.copy()
    
        pmax[:pmin.size] = pmax[:pmin.size] + pmin
    
        return Poly(pmax)

    __radd__ = __add__

    
    def __sub__(self, other):
        """Retorna la resta de dos polinomios. Si `other` es un número, se trata como un polinomio de grado cero `Poly([other])`"""
        if isinstance(other, (int, float, complex)):
            other = Poly([other])

        if self.deg < other.deg:
            pmin, pmax = self.coef, -other.coef.copy()
        else:
            pmin, pmax = -other.coef, self.coef.copy()
    
        pmax[:pmin.size] += pmin
    
        return Poly(pmax)

    
    def __rsub__(self, other):
        """Retorna la resta `other-self` de dos polinomios. Si `other` es un número, se trata como un polinomio de grado cero `Poly([other])`"""
        if isinstance(other, (int, float, complex)):
            other = Poly([other])

        if self.deg < other.deg:
            pmin, pmax = -self.coef, other.coef.copy()
        else:
            pmin, pmax = other.coef, -self.coef.copy()
    
        pmax[:pmin.size] += pmin
    
        return Poly(pmax)

    
    def __mul__(self, other):
        """Retorna la multiplicación entre dos polinomios. Si `other` es un número, se trata como un polinomio de grado cero `Poly([other])`"""
        if isinstance(other, (int, float, complex)):
            other = Poly([other])
        
        k = self.deg + 1
        coef = np.zeros(k + other.deg, dtype=self.coef.dtype)
        
        for i,c in enumerate(other.coef):
            coef[i:i+k] = coef[i:i+k] + c*self.coef
        
        return Poly(coef)

    __rmul__ = __mul__

    
    def deriv(self, m=1):
        """Retorna la m-ésima derivada del polinomio."""
        coef = [i*c for i,c in enumerate(self.coef[1:], 1)]

        # si coef es una lista vacia, retorna Poly([0])
        if not coef:
            coef = [0]

        p = Poly(coef)
        return p if (m==1 or coef==[0]) else p.deriv(m-1) 
    
    
    def __repr__(self):
        return f"PolyClass({self.coef})"


# esto es para testear el código sin usar import
if __name__ == "__main__":
    P = Poly([1,2,3,4])
    Q = Poly([1,1])

    print(f"P = {P}")
    print(f"Q = {Q}")
    print(f"2.5P = {2.5*P}")
    print(f"-P = {-P}")
    print(f"+P = {+P}")
    print(f"P+Q = {P+Q}")
    print(f"P-Q = {P-Q}")
    print(f"P(2.5) = {P(2.5)}")
    print(f"P([1, 2, 2.5]) = {P([1, 2, 2.5])}")

    print("\nDerivadas de P")
    for m in range(1, 6):
        print(f"derivada m={m}: {P.deriv(m)}")
