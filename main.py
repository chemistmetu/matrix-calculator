@ -0,0 +1,226 @@
import numpylike1 as mml
def find_all_eigenvalues(charac_eq):
    eigenvalues = []
    poly = charac_eq
    max_roots = len(poly.terms)
    attempts = max_roots * 4  
    for attempt in range(attempts):
        try:
            start = -max_roots + attempt * 0.5  
            root = mml.newton_raphson(poly, x0=start)
            root = round(root, 8)  
            if not any(abs(root - r) < 1e-6 for r in eigenvalues):
                eigenvalues.append(root)
                poly = mml.deflate(poly, root)
                if max(t.power for t in poly.terms) <= 1:
                    break
        except Exception:
            continue  

    max_power = max(t.power for t in poly.terms)
    if max_power == 2:
        coeffs = mml.coeffs_(poly)
        roots = mml.quadratic_formula(coeffs[0], coeffs[1], coeffs[2])
        if isinstance(roots, (list, tuple)):
            for r in roots:
                if not any(abs(r - ev) < 1e-6 for ev in eigenvalues):
                    eigenvalues.append(r)
        else:
            if not any(abs(roots - ev) < 1e-6 for ev in eigenvalues):
                eigenvalues.append(roots)
    elif max_power == 1:
        coeffs = mml.coeffs_(poly)
       
        root = -coeffs[1] / coeffs[0]
        if not any(abs(root - ev) < 1e-6 for ev in eigenvalues):
            eigenvalues.append(root)

    return eigenvalues

def deflate(poly, root, tol=1e-8):
    new_terms = []
    remainder = 0
    poly_terms_sorted = sorted(poly.terms, key=lambda t: -t.power)  

    n = poly_terms_sorted[0].power  
    new_coefs = [0]*(n)  

    
    coefs = [0]*(n+1)
    for t in poly_terms_sorted:
        coefs[n - t.power] = t.coef

    
    new_coefs[0] = coefs[0]
    for i in range(1, n):
        new_coefs[i] = coefs[i] + new_coefs[i-1]*root
    remainder = coefs[n] + new_coefs[n-1]*root

    if abs(remainder) > tol:
        raise ValueError(f"Deflation remainder not zero: {remainder}")

    
    new_terms = [mml.Term(c, n - 1 - i) for i, c in enumerate(new_coefs)]

    
    return mml.Polynomial(new_terms)









def menu():
    print("Matrix Calculator Menu")
    print("1. Construct Matrix")
    print("2. Addition of Matrices")
    print("3. Substraction of Matrices")
    print("4. Multiplication of Matrices")
    print("5. Multiplication of Single Matrix with a Scalar")
    print("6. Determinant of Matrix")
    print("7. Inverse process for Matrix")
    print("8. Transpose for the Matrix")
    print("9. Calculate eigenvalues and eigenvectors")
    print("0. Exit")

def print_matrix(M):
    for row in M:
        print(row)

def run():
    A = None
    B = None

    while True:
        menu()
        operation = input("Select Operation: ")

        if operation == "1":
            print("Constructing Matrix the first Matrix...")
            A = mml.const_matrix(int(input("Enter the Number of Rows for the first Matrix: ")), int(input("Enter the Number of Columns for the first Matrix: ")))
            B = None

        elif operation == "2":
            B = mml.const_matrix(int(input("Enter the Number of Rows for the second Matrix: ")), int(input("Enter the Number of Columns for the second Matrix:")))
            if A is None or B is None:
                print("Matrices are not defined.")
            else:
                try:
                    result = mml.add_matrix(A, B)
                    print("Result of Addition Process of the Matrices: ")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)

        elif operation == "3":
            B = mml.const_matrix(int(input("Enter the Number of Rows for the second Matrix: ")), int(input("Enter the Number of Columns for the second Matrix:")))
            if A is None or B is None:
                print("Matrices are not defined.")
            else:
                try:
                    result = mml.subs_matrix(A, B)
                    print("Result of Substraction Process: ")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)

        elif operation == "4":
            B = mml.const_matrix(int(input("Enter the Number of Rows for the second Matrix: ")), int(input("Enter the Number of Columns for the second Matrix: ")))
            if A is None or B is None:
                print("Matrices are not defined.")
            else: 
                try:
                    result = mml.mult_matrix(A, B)
                    print("Result of Multiplication of the Matrices: ")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)
      
        elif operation == "5":
            if A is None:
                print("Matrix is not defined.")
            else:
                try:
                    scalar = float(input("Enter value for scalar multiplication."))
                    result = mml.scalar_mult(A, scalar)
                    print(f"Result of Scalar Multiplication of the Matrix with {scalar}:")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)
        elif operation == "6":
            if A is None:
                print("Matrix is not defined.")
            else:
                try:
                    result = mml.det_matrix(A)
                    print(f"Determinant of the Matrix : {result}")
                except Exception as e:
                    print("Error:", e)

        elif operation == "7":
            if A is None:
                print("Matrix is not defined.")
            else:
                try:
                    result = mml.inverse_matrix(A)
                    print("Result of Inverse Process of the Matrix: ")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)

        elif operation == "8":
            if A is None:
                print("Matrix is not defined.")
            else:
                try:
                    result = mml.transpose(A)
                    print("Result of transpose Process of the Matrix: ")
                    print_matrix(result)
                except Exception as e:
                    print("Error:", e)
              





        elif operation == "9":
            if A is None:
              print("Matrix is not defined.")
            else:
                try:
                    x = mml.Polynomial([mml.Term(1, 1)])
                    charac_eq = mml.det_matrix(mml.A_minus_xI(A, x))
                    print("Characteristic Polynomial: ", charac_eq)
                          
                    eigenvalues = find_all_eigenvalues(charac_eq)
                    eigvecs = mml.eigen_vectors(A, eigenvalues)

                    print("Eigenvalues:")
                    for ev in eigenvalues:
                        print(ev)
                    

                    print("\nEigenvectors:")
                    for ev, vectors in eigvecs.items():
                        print(f"Eigenvalue: {ev}")
                        if vectors:
                            for v in vectors:
                                print(v)
                        else:
                            print("No eigenvectors is found.")

                except Exception as e:
                    print("Error:", e)
  
        elif operation == "0":
            print("Exiting from the program...")
            break
        else:
             print("Invalid try.")
     
if __name__ == "__main__":
    run()
