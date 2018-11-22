using Plots
gr(size=(600,600))
default(fmt=:png)

function declupivot(A::Matrix; diagtol = 1e-12)
    n = size(A, 2)
    p = collect(1:n)
    for j = 1:n-1
        # Quem é o pivô
        pivo, k = abs(A[j,j]), j
        for i = j+1:n
            if abs(A[i,j]) > pivo
                pivo, k = abs(A[i,j]), i
            end
        end
        if pivo <= diagtol
            error("Matriz singular ou muito próxima de ser singular")
        end

        if k != j
            p[k], p[j] = p[j], p[k]
            A[[k;j],:] = A[[j;k],:]
        end
        ajj = A[j,j]
        for i = j+1:n
            mij = A[i,j] / ajj
            A[i,j] = mij
            A[i,j+1:n] -= mij * A[j,j+1:n]
        end
    end
    L = tril(A, -1) + I
    U = triu(A)
    return p, L, U
end

function declurefine(A, b)
    p, L, U = declupivot(copy(A))
    c = b[p]
    y = L \ c
    x = U \ y
    E = zeros(length(b))
    println(p)
    println(x)
    println(A*x - b)
    for i = 1:length(b) #ver aqui
        xi = BigFloat.(x)
        ri = b - A * xi
        yi = L \ ri[p]
        Δi = U \ yi
        x = xi + Δi
        ei = BigFloat(norm(x - xi))
        E[i] = ei
        println("r$i = $(norm(ri))")
        println("x = $x")
        println("Erro = $ei")
    end
    return E
end

#A = [4.0 -1 0 -1; 1 -2 1 0; 0 4 -4 1; 5 0 5 -10]
#b = [6.0, 8, -7, -40]
#A = [2.0 3 -1; -3 5 6; 1 1 2]
#b = [-4.0, 19, 11]
#A = [1.0 1; 99.4 99.9]
#b = [1.0, 99.2]
#A = [2.0 1 7 4 -3 -1 4 4 7 0; 4 2 2 3 -2 0 3 3 4 1; 3 4 4 2 1 -2 2 1 9 -3; 9 3 5 1 0 5 6 -5 -3 4; 2 0 7 0 -5 7 1 0 1 6; 1 9 8 0 3 9 9 0 0 5; 4 1 9 0 4 3 7 -4 1 3; 6 3 1 1 6 8 3 3 0 2; 6 5 0 -7 7 -7 6 2 -6 1; 1 6 3 4 8 3 -5 0 -6 0]
#b = [86.0, 45, 52.5, 108, 66.5, 90.5, 139, 61, -43.5, 31]
#A = [4 -1 0 -1 0 0 0 0 0 0; -1 4 -1 0 -1 0 0 0 0 0; 0 -1 4 0 0 -1 0 0 0 0;
    #-1 0 0 4 -1 0 0 0 0 0; 0 -1 0 -1 4 -1 -1 0 0 0; 0 0 -1 0 -1 4 0 -1 0 0;
    #0 0 0 0 -1 0 4 -1 0 0; 0 0 0 0 0 -1 -1 4 -1 0; 0 0 0 0 0 0 0 -1 4 -1;
    #0 0 0 0 0 0 0 0 -1 4]
#b = [-110, -30, -40, -110, 0, -15, -90, -25, -55, -65]
E = declurefine(A,b)
D = E[1:10]
A=[1:10]
scatter(A, D, c=:blue, leg=:false, yaxis=:log)
name = @sprintf("Erro_refinamento4")
png(name)
