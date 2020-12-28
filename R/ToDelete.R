p = 8
x = rep(1, p)
y = rep(2, p)
Z = read.csv("C:\\Users\\armop\\Dropbox\\PHD\\Research\\NPNDAG\\Code\\Zmat.csv", header = TRUE, sep = ",")
S = read.csv("C:\\Users\\armop\\Dropbox\\PHD\\Research\\NPNDAG\\Code\\Smat.csv", header = TRUE, sep = ",")
L = read.csv("C:\\Users\\armop\\Dropbox\\PHD\\Research\\NPNDAG\\Code\\Lmat.csv", header = TRUE, sep = ",")


Z = as.matrix(Z)
S = as.matrix(S)
L = as.matrix(L)
P_0 = matrix(1/ p, p, p)
D = proj.ds(P_0, x, y, Z)$P

perm.proj(D = D, S = S,L = L, mu = 0.05,P_old = diag(1,p))
perm.proj()
