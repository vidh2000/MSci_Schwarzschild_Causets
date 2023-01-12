import numpy as np
import pandas as pd

def mymean (xs):
    return sum(xs)/len(xs)

def mystd (xs):
    xs2 = [xi*xi for xi in xs]
    N = len(xs)
    return np.sqrt(sum(xs2)/(N-1) - (sum(xs))**2/(N*(N-1)))

def std_combine(Ns, mus, stds):
    coeffs = []
    M = len(Ns)
    N = sum(Ns)
    for i in range(M):
        Ni = Ns[i]
        mui = mus[i]
        stdi = stds[i]
        term_i = (Ni-1)/(N-1)*stdi**2
        coeffs.append(term_i)
        for j in range(i, M):
            Nj = Ns[j]
            muj = mus[j]
            stdj = stds[j]
            term_mixed = Ni*Nj/(N*(N-1)) * (mui - muj)**2
            coeffs.append(term_mixed)
    return np.sqrt(sum(coeffs))


xs = [1,2,3,4,5]
ys = [2,3,4,5,6,7,8]
zs = [0.1, 19, 0.2, 7, 38, 39]


nx = len(xs)
ny = len(ys)
nz = len(zs)


X = mymean(xs)
Y = mymean(ys)
Z = mymean(zs)
print(f"\nThe means are:\n  -X = {X};\n  -Y = {Y};\n  -Z = {Z}.\n")


sx = round(mystd(xs), 7); npsx = round(np.std(xs, ddof = 1), 7)
sy = round(mystd(ys), 7); npsy = round(np.std(ys, ddof = 1), 7)
sz = round(mystd(zs), 7); npsz = round(np.std(zs, ddof = 1), 7)
print("\nThe stds are:")
print(f"  -sx: {sx}(mine),  {npsx}(numpy);")
print(f"  -sy: {sy}(mine),  {npsy}(numpy);")
print(f"  -sz: {sz}(mine),  {npsz}(numpy).")
print("(I was just checking I wrote the code correctly)")


sxy = round(mystd(xs+ys), 7)
sxz = round(mystd(xs+zs), 7)
syz = round(mystd(ys+zs), 7)
sxyz = round(mystd(xs+ys+zs), 7)

sxyc = round(std_combine([nx, ny], [X, Y], [sx, sy]), 7)
sxzc = round(std_combine([nx, nz], [X, Z], [sx, sz]), 7)
syzc = round(std_combine([ny, nz], [Y, Z], [sy, sz]), 7)
sxyzc = round(std_combine([nx, ny, nz], [X, Y, Z], [sx, sy, sz]), 7)


stds = pd.DataFrame([[sxy, sxyc],[sxz, sxzc], [syz, syzc], [sxyz, sxyzc]],
                  columns = ["correct, ", "combined"])
print("\nFINALLY TESTING STD COMBINATION\n", stds)






