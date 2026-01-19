@ -0,0 +1,114 @@
class Term:

    def __init__(self, coef, power):
        self.coef = coef
        self.power = power

    def __mul__(self, other):
        if isinstance(other, Term):
            return Term(self.coef * other.coef, self.power + other.power)
        elif isinstance(other, (int, float)):
            return Term(self.coef * other, self.power)
        else:
            raise NotImplementedError("Invalid Try.")

    def __rmul__(self, other):
        return self.__mul__(other)
 
    def __add__(self, other):
        if isinstance(other, Term) and self.power == other.power:
            return Term(self.coef + other.coef, self.power)
        else:
            raise NotImplementedError("Power of the terms should be same")

    def __neg__(self):
        return Term(-self.coef, self.power)

    def derivative(self):
        if self.power == 0:
            return Term(0,0)
        return Term(self.coef * self.power, self.power-1)

    def __repr__(self):
        if self.power == 0:
            return f"{self.coef}"
        elif self.power == 1:
            return f"{self.coef}*x"
        else:
            return f"{self.coef}*x^{self.power}"


class Polynomial:
    def __init__(self, terms = None):
        if terms is None:
            self.terms = []
        elif isinstance(terms, (int, float)):
            self.terms = [Term(terms, 0)]
        elif isinstance(terms, list):
            self.terms = terms
        else:
            raise TypeError("Invalid Try.")

    def add_term(self, term):
        for i, t in enumerate(self.terms):
            if t.power == term.power:
                self.terms[i] = t + term
                return
        self.terms.append(term)

    def __add__(self, other):
        if not isinstance(other, Polynomial):
            raise NotImplementedError
        result = Polynomial()
        for t in self.terms:
            result.add_term(t)
        for t in other.terms:
            result.add_term(t)
        return result

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            new_terms = [Term(t.coef * other, t.power) for t in self.terms]
            return Polynomial(new_terms)
        elif isinstance(other, Polynomial):
            result = Polynomial()
            for t1 in self.terms:
                for t2 in other.terms:
                    term_product = t1 * t2
                    result.add_term(term_product)
            return result
        else:
            raise NotImplementedError("Invalid Try.")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __sub__(self, other):
       if not isinstance(other, Polynomial):
           raise NotImplementedError("Invalid Try.")
       result = Polynomial()
       for t in self.terms:
           result.add_term(t)
       for t in other.terms:
           result.add_term(-t)
       return result

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            new_terms = [Term(t.coef / other, t.power) for t in self.terms]
            return Polynomial(new_terms)
        else:
            raise NotImplementedError("Invalid Try")

    def __repr__(self):
        return "+".join([str(t) for t in sorted(self.terms, key = lambda x: -x.power)])

    def derivative(self):
        deriv_terms = [t.derivative() for t in self.terms if t.coef != 0]
        return Polynomial(deriv_terms)

    def evaluate(self, val):
        result = 0
        for term in self.terms:
            result += term.coef * (val** term.power)
        return result
