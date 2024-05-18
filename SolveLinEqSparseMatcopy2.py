from sys import getsizeof
import numpy as np
import math

class SolveLinearEquationsPCG:
    def __init__(self, data, col, ptr, coorden, b, num_iter, tol):
        self.data = data
        # print("len data:",len(self.data))
        self.col = col
        # print("len col:",len(self.col))
        self.ptr = ptr
        # print("len ptr:",len(self.ptr))
        self.coords = coorden
        # print("len coords:",len(self.coords))
        self.b = b
        self.num_iter = num_iter
        self.tol = tol

    def pcg(self):
        x = np.zeros(len(self.b))
        r = self.b - self.dot(x)
        p = r
        resk = r
        norm_r = [np.linalg.norm(r)]

        for i in range(self.num_iter):
            Ap = self.dot(p)
            alpha = np.dot(r.T, r) / np.dot(p.T, Ap)
            x += alpha * p
            r = resk - alpha * Ap
            p = r + ((np.dot(r.T, r) / np.dot(resk.T, resk)) * p)
            resk = r
            norm_r.append(np.linalg.norm(r))

        return x, norm_r

    def jacobi(self, omega, u, res):
        # print(len(self.coords))
        coords = np.array(self.coords)
        # print(len(coords))
        # D = np.array([self.data[self.ptr[i]] for i in range(len(self.ptr)-1)])
        coords_Diagonalvalues =  [(i,coord) for i, coord in enumerate(coords) if coord[0] == coord[1]]
        # print("len coords_Diagonalvalues: ",len(coords_Diagonalvalues))
        # print("Coords index:", coords_Diagonalvalues)
        Diagmatrix=np.zeros(2368)
        # print("len diagmatrix: ",len(Diagmatrix))
        data=np.array(self.data)
        for i in range(len(coords_Diagonalvalues)-1):
            Diagmatrix[coords_Diagonalvalues[i][1][1]]=data[coords_Diagonalvalues[i][0]]

        D=Diagmatrix
        # print("len D: ",len(D))
        # print("len D: ",len(D))
        resnorm = np.linalg.norm(res)
        resloc = res.copy()
        delu = np.zeros(len(u))
        delu = omega * resloc / D
        # print("delu:",delu)
        resloc -= self.dot(delu)
        u += delu
        # u=1
        return u, resloc

    def iterate_jacobi(self, solini):
        sol = np.array(solini, dtype='float64')
        res = self.b - self.dot(sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.jacobi(1, sol, res)
            # print("res:",res)
            resnorm.append(np.linalg.norm(res))
        # print(self.coords)

        return sol, resnorm

    def SSOR(self, omega, u, res):
        delu, resc = self.gauss_seidel(u, res)
        delu2, resc2 = self.gauss_seidel_bb(delu, resc)
        return delu2, resc2

    def iterative_SSOR(self, solini, omega):
        sol = solini
        res = self.b - self.dot(sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.SSOR(omega, sol, res)
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

    def gauss_seidel(self, u, res):
        n = len(res)
        delu = [0] * n
        resc = res.copy()
        delu[0] = res[0] / self.data[0]

        for i in range(1, n):
            start = self.ptr[i]
            end = self.ptr[i + 1]
            resc[i] -= np.dot(self.data[start:end], delu[self.col[start:end]])
            delu[i] += resc[i] / self.data[self.ptr[i + 1] - 1]

        resc = res - self.dot(delu)
        return u + delu, resc

    def gauss_seidel_bb(self, u, res):
        n = len(res)
        delu = [0] * n
        resc = res.copy()
        delu[n-1] = res[n-1] / self.data[self.ptr[n-1]]

        for i in range(n-2, -1, -1):
            start = self.ptr[i]
            end = self.ptr[i + 1]
            resc[i] -= np.dot(self.data[start:end], delu[self.col[start:end]])
            delu[i] += resc[i] / self.data[self.ptr[i]]

        resc = res - self.dot(delu)
        return u + delu, resc

    def iterative_conjugate_gradient_PCG(self, solini, method):
        sol = solini
        zero = np.zeros(len(sol))
        res = self.b - self.dot(sol)
        z, _ = method(1, zero, res)
        zk = z
        p = z
        resk = res
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            alpha = np.dot(res.T, z) / np.dot(p.T, self.dot(p))
            sol = sol + alpha * p
            res = resk - alpha * self.dot(p)
            z, _ = method(1, zero, res)
            beta = np.dot(res.T, z) / np.dot(resk.T, zk)
            zk = z
            p = z + beta * p
            resk = res
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

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
        for coord in self.coords:
            [row,col]=coord
            # print(col)
            result[row] += self.data[i] * vector[col]
            i+=1
        return result
        # result = np.zeros(len(vector)-1)
        # print(max(self.coords))
        # for i in range(len(self.data)-1):
        #     result[self.coords[i][0]-1] += self.data[i] * vector[self.coords[i][1]-1]
        # return result

# Usage example


class CSRMatrix:
    def __init__(self, file_path, a):
        self.file_path = file_path
        self.a = a
        self.data = []
        self.col = []
        self.ptr = []
        self.coord = []

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
                    print("hoola2")
                    ok=0
                res=nterms%2368
                division=int(nterms//2368)
                if res==0:
                    res=2368
                    division-=1
                # print("res:",res)
                # print("division:",division)
                # print("nterms",nterms)
                if element != 0.0:
                
                        # division+=1
                    self.data.append(element)
                    self.col.append(res-1)
                    self.ptr.append(division)
                    



               
                # linha+=1
        # self.ptr.append(len(self.data))

        # self.ptr[0]=0
        # self.ptr[1]-=1
        # # self.ptr[-1]+=1
        # print("len data:",len(self.data))
        # print("len col:",len(self.col))
        # print("len ptr:",len(self.ptr))
        print("TERMINOS ",nterms)
        print("Primeros 10 terminos ")
        print(self.data[:20])
        print(self.col[:20])
        print( self.ptr[:20])
        print("Ultimos 10 terminos ")
        print(self.data[-20:])
        print(self.col[-20:])
        print( self.ptr[-20:])
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
        print("len data:",len(self.data))
        print("len col:",len(self.col))
        print("len ptr:",len(self.ptr))
        coords=[]
        for i in range(len(self.col)):
            coords.append([self.ptr[i], self.col[i]])
        
        self.coord = coords
        print("ultimo valor",coords[-1])
        print("primer valor",coords[0])

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
listvec=processall_linesvec(file_path, a)
print("len vectores: ",len(listvec))
b=np.array(listvec,dtype='float64')
num_iter = 10
tol = 1e-6
data = matrix.get_data()
# print(len(data))
col = matrix.get_col()
# print(col)
# print(len(col))
ptr = matrix.get_ptr()
# print(len(ptr))
coordenadas = matrix.get_coord()
print(coordenadas[-1])
# print(coordenadas[-2])
print(coordenadas[0])
# print(coordenadas[-1::-10])
# print(coords)

solver = SolveLinearEquationsPCG(data, col, ptr, coordenadas, b, num_iter, tol)
# print(coordenadas)
print(len(coordenadas))
print(len(data))
print(len(col))
lendiag=0
for i in range(len(coordenadas)):
    if coordenadas[i][0]==coordenadas[i][1]:
        print([data[i],coordenadas[i]])
        lendiag+=1
print("lendiag:",lendiag)

print(solver.iterate_jacobi(np.zeros(len(listvec), dtype='float64')))


# print(solver.pcg())

# print(solver.iterative_conjugate_gradient_PCG(np.array(len(listvec)), solver.SSOR))
# print(solver.iterative_conjugate_gradient_PCG(np.array(len(listvec)), solver.jacobi))
