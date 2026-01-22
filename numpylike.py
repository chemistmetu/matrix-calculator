@ -0,0 +1,263 @@
from sympylike import Term, Polynomial

def const_matrix(row, column):
    print(f" for {row} * {column} size matrix, enter your values: ")
    M = []
    for i in range(row):
        r_matrix = []
        for j in range(column):
            while True:
                try:
                    value = float(input(f"Characters [{i+1}, {j+1}]: "))
                    r_matrix.append(value)
                    break
                except ValueError:
                    print("Error occured, enter numerical input")
        M.append(r_matrix)
    if len(M) == row and all(len(r) == column for r in M):
        print("\n The size of matrix is correct")
    else:
        print("\n The size of matrix is incorrect")
    return M

def add_matrix(A, B):
    if len(A) != len(B) or not all(len(rowA) == len(rowB) for rowA, rowB in zip(A, B)):
        raise ValueError("Error occured, check the lengths. ")
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def subs_matrix(A, B):
    if len(A) != len(B) or not all(len(rowA) == len(rowB) for rowA, rowB in zip(A, B)):
        raise ValueError("Error occured, check the lengths. ")
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))] for i in range(len(A))]

def scalar_mult(A, scalar):
    return [[A[i][j] * scalar for j in range(len(A[0]))] for i in range(len(A))] 

def transpose(A):
    rows = len(A)
    cols = len(A[0])
    return [[A[i][j] for i in range(rows)] for j in range(cols)]

def mult_matrix(A, B):
    if len(A[0]) != len(B):
        raise ValueError(" number of columns of A = number of rows of B.")
    result = []
    for i in range(len(A)):
        row = []
        for j in range(len(B[0])):
            overall = 0
            for k in range(len(B)):
                overall += A[i][k] * B[k][j]
            row.append(overall)
        result.append(row)
    return result

def det_matrix(A):
    if len(A) != len(A[0]):
        raise ValueError("Determinant process require square matrix.")
    def poly(x):
        if isinstance(x, Polynomial):
            return x
        return Polynomial([Term(x, 0)])
    n = len(A)
    if n == 1:
        return poly(A[0][0])
    if n == 2:
        return poly(A[0][0]) * poly(A[1][1]) - poly(A[0][1]) * poly(A[1][0]) 
    det = Polynomial(0)
    for c in range(n):
        sub = [row[:c] + row[c+1:] for row in A[1:]]
        cofactor = ((-1)**c) * poly(A[0][c]) * det_matrix(sub)
        det += cofactor
    return det

def poly_to_scalar(p):
    if len(p.terms) == 1 and p.terms[0].power == 0:
        return p.terms[0].coef
    else:
        raise ValueError("Invalid type")

def inverse_matrix(A):
    det = det_matrix(A)
    try:
        det_scalar = poly_to_scalar(det)
    except ValueError:
        raise ValueError("Invalid Try.")
    if det_scalar == 0:
        raise ValueError("This specific matrix has no inverse property.")
    n = len(A)
    if n == 2:
        return [[A[1][1]/det_scalar, -A[0][1]/det_scalar],
                [-A[1][0]/det_scalar, A[0][0]/det_scalar]]
    cofactors = []
    for r in range(n):
        cofactor_row = []
        for c in range(n):
            minor = [row[:c] + row[c+1:] for row in (A[:r] + A[r+1:])]
            cofactor_row.append(((-1)**(r+c)) * det_matrix(minor))
        cofactors.append(cofactor_row)
    cofactors = transpose(cofactors)
    for r in range(n):
        for c in range(n):
            cofactors[r][c] = cofactors[r][c]/det_scalar
    return cofactors

def my_eye_like(rowA, columnA):
    if rowA != columnA:
        raise ValueError("Invalid Try.")
    return [[1.0 if i == j else 0.0 for j in range(columnA)] for i in range(columnA)]

def A_minus_xI(A, x):
    n = len(A)
    I = my_eye_like(n, n)
    M = []
    for i in range(n):
        row = []
        for j in range(n):
            a_ij = Polynomial([Term(A[i][j], 0)])
            if i == j:
                val = a_ij - x
            else:
                val = a_ij 
            row.append(val)
        M.append(row)
    return M 

def newton_raphson(poly, x0 = 0.0, tol = 1e-10, max_iter = 1000):
    x = x0
    for i in range(max_iter):
        fx = poly.evaluate(x)
        f_pr = poly.derivative().evaluate(x)
        if abs(f_pr) < 1e-14:
            raise ValueError("derivative of the function is very small.")
        x_new = x - fx / f_pr
        if abs(x_new - x) < tol:
            return x_new
        x = x_new
    raise ValueError("The iterations failed.")

def deflate(poly, root):
    terms = sorted(poly.terms, key = lambda t: t.power, reverse = True)
    new_terms = []
    remainder = 0
    for i, term in enumerate(terms):
        if i == 0:
            coef = term.coef
            new_terms.append(Term(coef, term.power - 1))
        else:
            coef = term.coef + remainder * root
            if term.power > 0:
                new_terms.append(Term(coef, term.power - 1))
            remainder = coef
    return Polynomial(new_terms)

def coeffs_(poly):
    max_power = max(term.power for term in poly.terms)
    coefficients = [0] * (max_power + 1)
    for term in poly.terms:
        coefficients[term.power] = term.coef
    coefficients.reverse()
    return coefficients

def quadratic_formula(a, b, c):
    discriminant = b**2 - 4*a*c
    if discriminant > 0:
        root1 = (-b + discriminant**0.5) / (2*a)
        root2 = (-b - discriminant**0.5) / (2*a)
        return root1, root2
    elif discriminant == 0:
        root_d = -b / (2*a)
        return (root_d,)
    else:
        real = -b / (2*a)
        imag = (abs(discriminant)**0.5) / (2*a)
        return complex(real, imag), complex(real, -imag)

def rref(matrix, tol=1e-12):
    
    A = [ [complex(x) for x in row] for row in matrix ]  
    rows = len(A)
    cols = len(A[0]) if rows>0 else 0
    r = 0
    for c in range(cols):
        if r >= rows:
            break
       
        pivot = None
        max_val = 0
        for i in range(r, rows):
            val = abs(A[i][c])
            if val > max_val and val > tol:
                max_val = val
                pivot = i
        if pivot is None:
            continue
     
        A[r], A[pivot] = A[pivot], A[r]
        pivot_val = A[r][c]
        
        A[r] = [elem / pivot_val for elem in A[r]]
      
        for i in range(rows):
            if i != r:
                factor = A[i][c]
                if abs(factor) > tol:
                    A[i] = [A[i][j] - factor * A[r][j] for j in range(cols)]
        r += 1
   
    for i in range(rows):
        for j in range(cols):
            if abs(A[i][j]) < tol:
                A[i][j] = 0+0j
    return A

def null_space(matrix, tol=1e-9):
    
    R = rref(matrix, tol=tol)
    rows = len(R)
    cols = len(R[0]) if rows>0 else 0

    pivot_cols = []
    for i in range(rows):
       
        pivot_col = None
        for j in range(cols):
            if abs(R[i][j]) > tol:
                pivot_col = j
                break
        if pivot_col is not None:
            pivot_cols.append(pivot_col)

    free_vars = [j for j in range(cols) if j not in pivot_cols]
    if not free_vars:
        return [] 

    basis = []
    for free in free_vars:
        vec = [0+0j] * cols
        vec[free] = 1+0j
        
        for i in range(rows):
           
            pivot_col = None
            for j in range(cols):
                if abs(R[i][j] - 1) < tol:   
                    pivot_col = j
                    break
            if pivot_col is not None:
               
                vec[pivot_col] = -(R[i][free])
       
        basis.append([ (v.real if abs(v.imag) < tol else v) for v in vec ])
    return basis


def eigen_vectors(A, eigenvalues):
    n = len(A)
    I = my_eye_like(n, n)
    eigvecs = {}
    for ev in eigenvalues:
        M = [[A[i][j] - ev * I[i][j] for j in range(n)] for i in range(n)]
        ns = null_space(M)
        eigvecs[ev] = ns
    return eigvecs
