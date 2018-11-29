using Plots
gr(size=(600,600))
default(fmt=:png)

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
    E = zeros(length(b)*2)
    println(p)
    println(x)
    println(A*x - b)
    for i = 1:length(b)*2 
        xᵢ = BigFloat.(x)
        rᵢ = b - A * xᵢ
        yᵢ = L \ rᵢ[p]
        Δᵢ = U \ yᵢ
        x = xᵢ + Δᵢ
        eᵢ = BigFloat(norm(x - xᵢ))
        E[i] = eᵢ
        println("r$i = $(norm(rᵢ))")
        println("x = $x")
        println("Erro = $eᵢ")
    end
    return E
end

E = declurefine(A,b)
D = E[1:10]
A=[1:10]
scatter(A, D, c=:blue, leg=:false, yaxis=:log)
name = @sprintf("Errorefinamento")
png(name)
