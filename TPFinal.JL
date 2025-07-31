using FASTX
using Plots

# -------------------------------
# Función para cargar una secuencia desde un archivo FASTA
function cargar_secuencia_fasta(path::String)::String
    open(FASTA.Reader, path) do reader
        for record in reader
            return String(sequence(record))
        end
    end
end

# -------------------------------
# Alineamiento global - Needleman-Wunsch
function needleman_wunsch(seq1, seq2; match=1, mismatch=-1, gap=-1)
    m, n = length(seq1), length(seq2)
    score = fill(0, m+1, n+1)

    for i in 1:m+1
        score[i, 1] = (i-1) * gap
    end
    for j in 1:n+1
        score[1, j] = (j-1) * gap
    end

    for i in 2:m+1, j in 2:n+1
        diag = score[i-1, j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)
        up = score[i-1, j] + gap
        left = score[i, j-1] + gap
        score[i, j] = maximum([diag, up, left])
    end

    return score[m+1, n+1], score
end

# -------------------------------
# Alineamiento local - Smith-Waterman
function smith_waterman(seq1, seq2; match=2, mismatch=-1, gap=-2)
    m, n = length(seq1), length(seq2)
    score = fill(0, m+1, n+1)
    max_score = 0

    for i in 2:m+1, j in 2:n+1
        diag = score[i-1, j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)
        up = score[i-1, j] + gap
        left = score[i, j-1] + gap
        score[i, j] = max(0, diag, up, left)
        max_score = max(max_score, score[i, j])
    end

    return max_score, score
end

# -------------------------------
# Dotplot simple por subsecuencias comunes
# Función para generar el dotplot entre dos secuencias
function dotplot(seq1::String, seq2::String; window::Int=4)
    dict1 = Dict{String, Vector{Int}}()
    dict2 = Dict{String, Vector{Int}}()

    for i in 1:(length(seq1) - window + 1)
        subseq = seq1[i:i+window-1]
        push!(get!(dict1, subseq, Int[]), i)
    end

    for j in 1:(length(seq2) - window + 1)
        subseq = seq2[j:j+window-1]
        push!(get!(dict2, subseq, Int[]), j)
    end

    # Encuentra coincidencias exactas de subsecuencias
    matches = intersect(keys(dict1), keys(dict2))

    x = Int[]
    y = Int[]

    for match in matches
        for i in dict1[match], j in dict2[match]
            push!(x, i)
            push!(y, j)
        end
    end

    # Graficar los puntos
    scatter(x, y,
        xlabel = "Posición en Secuencia 1",
        ylabel = "Posición en Secuencia 2",
        title = "Dotplot (window = $window)",
        markersize = 2,
        legend = false,
        grid = true)
end

# -------------------------------
# MAIN - Ejecución
seq1 = cargar_secuencia_fasta("./Fasta/1VXA.fasta")
seq2 = cargar_secuencia_fasta("./Fasta/1MBN.fasta")

puntaje_global, _ = needleman_wunsch(seq1, seq2)
println("\nNeedleman-Wunsch (Global): Puntaje = ", puntaje_global)

puntaje_local, _ = smith_waterman(seq1, seq2)
println("Smith-Waterman (Local): Puntaje = ", puntaje_local)

# Generar dotplot con ventana de 4 (o cambiar a gusto)
window = 3
dotplot(seq1, seq2, window=window)