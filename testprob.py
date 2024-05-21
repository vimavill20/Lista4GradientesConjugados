import numpy as np
class SolveLinearEquationsPCG:
    def __init__(self, data, col, ptr, coorden, b, num_iter, tol):
        self.data = data
        print("len data:",len(self.data))
        self.col = col
        print("len col:",len(self.col))
        self.ptr = ptr
        print("len ptr:",len(self.ptr))
        self.coords = coorden
        print("len coords:",len(self.coords))
        self.b = b
        print("len b:",len(self.b))
        self.num_iter = num_iter
        self.tol = tol

    def pcg(self,solini):
        # x = np.zeros(len(self.b))
        x=solini.copy()
        r = self.b - self.dot(x)
        # print(np.sum(r)==np.sum(x))
        p = r
        resk = r
        norm_r = [np.linalg.norm(r)]

        for i in range(self.num_iter):
            Ap = self.dot(p)
            alpha = np.dot(r.T, r) / np.dot( self.dot(p),p.T)
            x += alpha * p
            r = resk - (alpha * self.dot(p))

            p = r + ((np.dot(r.T, r) / np.dot(resk.T, resk)) )
            resk = r.copy()
            norm_r.append(np.linalg.norm(r))

        return x, norm_r

    def jacobi(self, omega, u, res):
        # print("HOLA")
        D = self.get_diag_matrix()
        resnorm = np.linalg.norm(res)
        resloc = res.copy()
        delu = np.zeros(len(u))
        delu =  omega *np.divide(resloc, D)
        resloc -= self.dot(delu)
        
        return u+delu, resloc

    def iterate_jacobi(self, solini):
        sol = np.array(solini, dtype='float64')
        res = self.b - self.dot(sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.jacobi(1.0, sol, res)
            resnorm.append(np.linalg.norm(res))
            # if resnorm[-1] >= resnorm[-2]:  # Changed > to >= for stability
            #     break

        return sol, resnorm

    

    def gauss_seidel(self, u, res):
        n = len(res)
        delu = np.zeros(n)
        resc = res.copy()
        omega=1
        delu[0] = omega * res[0] / self.Aij(0, 0)

        for i in range(1, n):
            lower_sum = self.InnerProductLowerRows(delu, i)
            resc[i] -= lower_sum
            delu[i] += omega * resc[i] / self.Aij(i, i)

        resc = res - self.dot(delu)
        return u + delu, resc

    def gauss_seidel_bb(self, u, res):
        
        n = len(res)
        omega=1
        delu = np.zeros(n)
        resc = res.copy()

        delu[n - 1] = omega * res[n - 1] / self.Aij(n - 1, n - 1)

        for i in range(n - 2, -1, -1):
            upper_sum = self.InnerProductUpperRows(delu, i)
            resc[i] -= upper_sum
            delu[i] += omega * resc[i] / self.Aij(i, i)

        resc = res - self.dot(delu)
        return u + delu, resc
    def SSOR(self, omega, u, res):
        delu, res = self.gauss_seidel(u, res)
        delu, res = self.gauss_seidel_bb(delu, res)
        return delu, res

    def iterative_SSOR(self, solini, omega):
        sol = solini
        res = self.b - self.dot(sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.SSOR(omega, sol, res)
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

    def iterative_conjugate_gradient_PCG(self, solini, method):
        sol = solini
        zero = np.zeros(len(sol))
        sol = zero
        res = self.b - self.dot(sol)
        z, _ = method(1, zero, res)
        zk = z
        p = z
        resk = res
        resnorm = [np.linalg.norm(res)]
        print("iteraciones",self.num_iter)
        for i in range(self.num_iter):
            Ap= self.dot(p)
            alpha = np.dot(res, z) / np.dot(p, self.dot(p))
            sol = sol + alpha * p
            res = resk - alpha * Ap
            z, _ = method(1, zero, res)
            beta = np.dot(res, z) / np.dot(resk, zk)
            zk = z
            p = z + beta * p
            resk = res
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm
    def Aij(self, i, j):
        # if any ([i,j] == coords for coords in self.coords):
        #     print("HOLA1")

        #     return self.data[j]
        # print("HOL2")
        for k in range(self.ptr[i], self.ptr[i+1]):
            if self.col[k] == j:
                return self.data[k]
        
        return 0.0

    # def dot(self, vector):
    #     if len(vector) != len(self.ptr)-1:
    #         raise ValueError("Incompatible dimensions for matrix-vector multiplication")

    #     result = [0] * (len(self.ptr)-1)
    #     for i in range(len(self.ptr) - 1):
    #         start = self.ptr[i]
    #         end = self.ptr[i + 1]
    #         for j in range(start, end):
    #             result[i] += self.data[j] * vector[self.col[j]]
    #             # Use coords array to access the column index
    #             col_index = self.coords[self.col[j]][1]
    #             result[i] += self.data[j] * vector[col_index]
    #     return result
    # def dot(self, vector):
    #     result = [0] * len(coords)
    #     for i in range(len(ptr) - 1):
    #         for j in range(ptr[i], ptr[i+1]):
    #             result[i] += data[j] * vector[coords[col[j]]]
    #     return result
    def dot(self, vector):
        result = np.zeros(len(vector))
        i=0
        # print("HOLA V")
        for coord in self.coords:
            # print(len(self.coords))
            # print(len(coord))
            [row,col]=coord
            # print(col)
            result[row] += self.data[i] * vector[col]
            i+=1
        return result
    def InnerProductLowerRows(self, vector: np.array, row: int) -> float:
        """
        Multiply the lower part of the sparse matrix by a vector.
        """
        nelem = self.ptr[row + 1] - self.ptr[row]
        ntotal = self.ptr[row]

        sum = 0
        for j in range(ntotal, ntotal + nelem):
            if self.col[j] >= row:
                continue

            sum += self.data[j] * vector[self.col[j]]

        return sum

    def InnerProductUpperRows(self, vector: np.array, row: int) -> float:
        """
        Multiply the upper part of the sparse matrix by a vector.
        """
        nelem = self.ptr[row + 1] - self.ptr[row]
        ntotal = self.ptr[row]

        sum = 0
        for j in range(ntotal, ntotal + nelem):
            if self.col[j] <= row:
                continue

            sum += self.data[j] * vector[self.col[j]]

        return sum
    # def dot(self, vector):
    #     print("HOLA C")
    #     prod = np.zeros(len(self.ptr) - 1)
    #     for i in range(len(self.ptr) - 1):
    #         sum = 0.0
    #         for k in range(self.ptr[i], self.ptr[i+1]):
    #             sum += self.data[k] * vector[self.col[k]]

    #         prod[i] = sum

    #     return prod
        
    def get_diag_matrix(self):
        DiagMatrix = np.zeros(len(self.ptr)-1)
        indexDiagMatrix = 0
        for i in range(len(self.coords)):
            if self.coords[i][0] == self.coords[i][1]:
                DiagMatrix[indexDiagMatrix] = self.data[i]
                indexDiagMatrix += 1
        return DiagMatrix
class CSRMatrix:
    def __init__(self, file_path, a):
        self.file_path = file_path
        self.a = a
        self.data = []
        self.col = []
        self.ptr = []
        self.coord = []
        self.row = []
        # self.dim=dim

    def processall_lines(self):
        with open(self.file_path, 'r') as file:
            lines = file.readlines()
        
        self.data = []
        self.col = []
        self.ptr = []
        col = 0
        current_ptr = 0
        linha=0
        col=[]
        row=[]
        rowindex=0
        nterms=0
        fake=0
        for line in lines[self.a:]:
            line = line.replace('{', '')
            line = line.replace('}', '')
            elements = [float(element) for element in line.split(',') if element.strip()]
            for element in elements:
                nterms+=1
                fake+=1
                if nterms==2366:
                    # print("hoola2")
                    ok=0
                res=nterms%2368
                division=int(nterms//2368)
                if res==0:
                    res=2368
                    division-=1
                    self.ptr.append(len(self.data))
                # print("res:",res)
                # print("division:",division)
                # print("nterms",nterms)
                if element != 0.0:
                        # division+=1
                    self.data.append(element)
                    self.col.append(res-1)
                    self.row.append(division)
                    
        self.ptr.insert(0,0)
        self.ptr[-1]+=1
        # print("SUMA DE ELEMENTOS",sum(self.ptr))


               
                # linha+=1
        # self.ptr.append(len(self.data))

        # self.ptr[0]=0
        # self.ptr[1]-=1
        # # self.ptr[-1]+=1
        # print("len data:",len(self.data))
        # print("len col:",len(self.col))
        # print("len ptr:",len(self.ptr))
        # print("TERMINOS ",nterms)
        # print("Primeros 10 terminos ")
        # print(self.data[:20])
        # print(self.col[:20])
        # print( self.ptr[:20])
        # print("Ultimos 10 terminos ")
        # print(self.data[-20:])
        # print(self.col[-20:])
        # print( self.ptr[-20:])
        # # self.ptr.append(len(self.data)+1)
        # nnzrow = np.diff(self.ptr)
        # print("nnzrow:",nnzrow)
        # print(len(nnzrow))
        # print(nnzrow[:10])
        # print("suma",sum(nnzrow))
        # coords = []
        # rows = np.repeat(np.arange(len(self.ptr)-1), np.diff(self.ptr))
        # print(rows[:17])
        # print(len(rows))
        #######             ########
        #####COMPROBACIONES ########
        #####               ########
        ########            ########
        # print("len data:",len(self.data))
        # print("len col:",len(self.col))
        # print("len ptr:",len(self.ptr))
        coords=[]
        for i in range(len(self.col)):
            coords.append([self.row[i], self.col[i]])
        
        self.coord = coords
        # print("ultimo valor",coords[-1])
        # print("primer valor",coords[0])
    
    def get_data(self):
        return self.data

    def get_col(self):
        return self.col

    def get_ptr(self):
        return self.ptr

    def get_coord(self):
        return self.coord
    
    # def dot(self, vector):
    #     if len(vector) != len(self.matrix[0]):
    #         raise ValueError("Dimensiones incompatibles para multiplicacion matriz-vector")

    #     result = [0] * len(self.matrix)
    #     for i in range(len(self.matrix)):
    #         start = self.ptr[i]
    #         end = self.ptr[i + 1]
    #         for j in range(start, end):
    #             result[i] += self.data[j] * vector[self.col[j]]
    #     return result
    def dot(self, vector):
        result = np.zeros(len(vector))
        i=0
        for coord in self.coords:
            [row,col]=coord
            result[row-1] += self.data[i] * vector[col]
            i+=1
        return result
    


# Usage example
matrix = CSRMatrix('/Users/victorvillegassalabarria/Downloads/matrix.txt', 1)
# matrix = CSRMatrix('/Users/victorvillegassalabarria/Downloads/matriztest.txt', 0)

matrix.processall_lines()

#2368x2368
def processall_linesvec(file_path, a):
        with open(file_path, 'r') as file:
            lines = file.readlines()
        b = len(lines)
        vec=[]
        for i in range(a, b):
            lines[i] = lines[i].replace('{', '')
            lines[i] = lines[i].replace('\n', '')
            lines[i] = lines[i].replace('}', '')
            a = lines[i].split(',')

            for element in a:
                if element != ' ':
                    if element != '':
                        element = float(element)
                        vec.append(element)
        return vec

file_path='/Users/victorvillegassalabarria/Downloads/rhs.txt'
a=1
import sys 
# sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

def processall_linesvec(file_path):
    with open(file_path, "r") as f:
        vector = np.array([float(x) for x in f.read().split(",")])      
    
    return vector

listvec=processall_linesvec(file_path)
print("len vectores: ",len(listvec))
b=np.array(listvec,dtype='float64')
num_iter = 500
tol = 1e-15
data = list(matrix.get_data())
# print(len(data))
col = np.array(matrix.get_col())
# print(col)
# print(len(col))
ptr = np.array(matrix.get_ptr())
# print(len(ptr))
coordenadas = list(matrix.get_coord())
DiagMatrix = np.zeros(len(ptr)-1)
indexDiagMatrix=0
lendiag=0
for i in range(len(coordenadas)):
    if coordenadas[i][0]==coordenadas[i][1]:
        if coordenadas[i][0]<=100:
            print([data[i],coordenadas[i]])
        DiagMatrix[indexDiagMatrix]=data[i]
        lendiag+=1
        indexDiagMatrix+=1

solver = SolveLinearEquationsPCG(data, col, ptr, coordenadas, b, num_iter, tol)
# jacsol=solver.iterate_jacobi(np.zeros(len(b)))
# jacsol=solver.pcg(np.zeros(len(b)))
# jacsol=solver.iterative_SSOR(np.zeros(len(b)), 1)
jacsol=solver.iterative_conjugate_gradient_PCG(np.zeros(len(b)), solver.SSOR)
# jacsol=solver.iterative_conjugate_gradient_PCG(np.zeros(len(b)), solver.jacobi)

print(jacsol[0])
print(jacsol[1])
