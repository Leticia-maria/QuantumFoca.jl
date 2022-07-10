function computeenergy(basis, molecule::Molecule, maxiter=20, convergence=1e-6)
    S = overlap(basis, molecule)
    T = kinetic(basis, molecule)
    V = attraction(basis, molecule)
    G = repulsion(basis, molecule)

    K = length(basis)

    Hcore = T .+ V

    D = zeros(K, K)
    P = zeros(K, K)

    X = sqrt(inv(S))

    Eel = 0.0

    for iteration in 0:maxiter
        Eold = Eel
        for n in 1:K
            for m in 1:K
                P[m, n] = 0.0
                for ℓ in 1:K
                    for s in 1:K
                        P[m, n] += D[ℓ, s] * (G[m, n, s, ℓ] - 0.5 * G[m, ℓ, s, n])
                    end
                end
            end
        end

        F = Hcore + P 
        Fp = X .* F .* X
        Cp = eigvecs(Fp)
        C = X .* Cp

        for n in 1:K
            for m in 1:K
                D[m, n] = 0.0
                for a in 1:trunc(Int64, N/2)
                    D[m, n] += 2 * (C[m, a] * C[n, a])
                end
            end
        end

        Eel = 0.0

        for m in 1:K
            for n in 1:K
                Eel += 0.5 * D[n, m] * (Hcore[m, n] + F[m, n])
            end
        end

        if (abs(Eel - Eold) < convergence) && (iteration > 0)
            break
        else
            continue
        end
    end