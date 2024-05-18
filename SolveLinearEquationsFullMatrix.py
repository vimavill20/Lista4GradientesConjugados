import numpy as np

class SolveLinearEquationsPCG:
    def __init__(self, A, b, num_iter, tol):
        self.A = A
        self.b = b
        self.num_iter = num_iter
        self.tol = tol

    def pcg(self):
        x = np.zeros(len(self.b))
        r = self.b - np.dot(self.A, x)
        p = r
        resk = r
        norm_r = [np.linalg.norm(r)]

        for i in range(self.num_iter):
            Ap = np.dot(self.A, p)
            alpha = np.dot(r.T, r) / np.dot(p.T, Ap)
            x += alpha * p
            r = resk - alpha * Ap
            p = r + ((np.dot(r.T, r) / np.dot(resk.T, resk)) * p)
            resk = r
            norm_r.append(np.linalg.norm(r))

        return x, norm_r

    def jacobi(self, omega, u, res):
        D = np.diag(self.A)
        resnorm = np.linalg.norm(res)
        resloc = res.copy()
        delu = np.zeros(len(u))
        delu = omega * resloc / D
        resloc -= np.dot(self.A, delu)
        u += delu

        return u, resloc

    def iterate_jacobi(self, solini):
        sol = solini
        res = self.b - np.dot(self.A, sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.jacobi(1, sol, res)
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

    def SSOR(self, omega, u, res):
        delu, resc = self.gauss_seidel(u, res)
        delu2, resc2 = self.gauss_seidel_bb(delu, resc)
        return delu2, resc2

    def iterative_SSOR(self, solini, omega):
        sol = solini
        res = self.b - np.dot(self.A, sol)
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            sol, res = self.SSOR(omega, sol, res)
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

    def gauss_seidel(self, u, res):
        n = len(res)
        delu = [0] * n
        resc = res.copy()
        delu[0] = res[0] / self.A[0][0]

        for i in range(1, n):
            resc[i] -= np.dot(self.A[i][0:i], delu[0:i])
            delu[i] += resc[i] / self.A[i][i]

        resc = res - np.dot(self.A, delu)
        return u + delu, resc

    def gauss_seidel_bb(self, u, res):
        n = len(res)
        delu = [0] * n
        resc = res.copy()
        delu[n-1] = res[n-1] / self.A[n-1, n-1]

        for i in range(n-2, -1, -1):
            resc[i] -= np.dot(self.A[i, i+1:], delu[i+1:])
            delu[i] += resc[i] / self.A[i, i]

        resc = res - np.dot(self.A, delu)
        return u + delu, resc

    def iterative_conjugate_gradient_PCG(self, solini, method):
        sol = solini
        zero = np.zeros(len(sol))
        res = self.b - np.dot(self.A, sol)
        z, _ = method(1, zero, res)
        zk = z
        p = z
        resk = res
        resnorm = [np.linalg.norm(res)]

        for i in range(self.num_iter):
            alpha = np.dot(res.T, z) / np.dot(p.T, np.dot(self.A, p))
            sol = sol + alpha * p
            res = resk - alpha * np.dot(self.A, p)
            z, _ = method(1, zero, res)
            beta = np.dot(res.T, z) / np.dot(resk.T, zk)
            zk = z
            p = z + beta * p
            resk = res
            resnorm.append(np.linalg.norm(res))

        return sol, resnorm

# Exemplo de uso
A = np.array([[4, 3, 0], [3, 4, -1], [0, -1, 4]])
b = np.array([24, 30, -24])
num_iter = 3
tol = 1e-6

solver = SolveLinearEquationsPCG(A, b, num_iter, tol)
print(solver.pcg())
print(solver.iterative_SSOR(np.array([0.0, 0.0, 0.0]), 1))
print(solver.iterative_conjugate_gradient_PCG(np.array([0.0, 0.0, 0.0]), solver.SSOR))
print(solver.iterative_conjugate_gradient_PCG(np.array([0.0, 0.0, 0.0]), solver.jacobi))

