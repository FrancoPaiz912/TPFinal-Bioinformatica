using FASTX
using BioSequences
using BioAlignments
using Plots
using BioSequences: LongAA

const aa_validos = Set(['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'])

function limpiar_secuencia(seq::AbstractString)
    return join([c for c in uppercase(seq) if c in aa_validos])
end

function cargar_fasta_completo(path::String)
    secuencias = Vector{Tuple{String, LongAA}}()
    open(FASTA.Reader, path) do reader
        for record in reader
            header = FASTX.description(record)
            raw_seq = String(FASTX.sequence(record))
            seq_limpia = limpiar_secuencia(raw_seq)
            push!(secuencias, (header, LongAA(seq_limpia)))
        end
    end
    return secuencias
end

function buscar_secuencia_por_texto(secuencias::Vector{Tuple{String, LongAA}}, texto::String)
    texto_lower = lowercase(texto)
    for (header, seq) in secuencias
        if occursin(texto_lower, lowercase(header))
            return (header, seq)
        end
    end
    error("No se encontró secuencia que contenga '$texto'")
end

function alinear_y_mostrar(seq1::LongAA, seq2::LongAA, nombre::String)
    modelo = AffineGapScoreModel(match=1, mismatch=-1, gap_open=-2, gap_extend=-1)

    aln_global = pairalign(GlobalAlignment(), seq1, seq2, modelo)
    println("\n===== Alineamiento global (Needleman-Wunsch) - $nombre =====")
    println(aln_global)

    aln_local = pairalign(LocalAlignment(), seq1, seq2, modelo)
    println("\n===== Alineamiento local (Smith-Waterman) - $nombre =====")
    println(aln_local)
end

# Función para extraer parte del header desde '>' hasta el primer '|'
function extraer_nombre(header::String)
    idx = findfirst('|', header)
    if idx === nothing
        return header
    else
        return header[1:idx-1]
    end
end

function dotplot(seq1::String, seq2::String, name1::String, name2::String; window::Int=4)
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

    matches = intersect(keys(dict1), keys(dict2))
    x = Int[]
    y = Int[]

    for match in matches
        for i in dict1[match], j in dict2[match]
            push!(x, i)
            push!(y, j)
        end
    end

    scatter(x, y,
        xlabel = "Posición en $name1",
        ylabel = "Posición en $name2",
        title = "Dotplot (window = $window)",
        markersize = 2,
        legend = false,
        grid = true,
        yflip = true)

        display(current())
end

# --- Comparaciones ---
comparaciones = [
    (file1 = "./Fasta/rcsb_pdb_3I40.fasta", texto1 = "Insulin A chain",
     file2 = "./Fasta/rcsb_pdb_2ZP6.fasta", texto2 = "Insulin A chain",
     nombre = "Insulina humana (Cadena A) vs bovina (Cadena A)"),

    (file1 = "./Fasta/rcsb_pdb_3I40.fasta", texto1 = "Insulin B chain",
     file2 = "./Fasta/rcsb_pdb_2ZP6.fasta", texto2 = "Insulin B chain",
     nombre = "Insulina humana (Cadena B) vs bovina (Cadena B)"),

    (file1 = "./Fasta/rcsb_pdb_1A3N.fasta", texto1 = "HEMOGLOBIN (ALPHA CHAIN)",
     file2 = "./Fasta/rcsb_pdb_2MBW.fasta", texto2 = "MYOGLOBIN",
     nombre = "Hemoglobina alfa humana vs Mioglobina cachalote"),

    (file1 = "./Fasta/rcsb_pdb_3EXJ.fasta", texto1 = "p53 core domain",
     file2 = "./Fasta/rcsb_pdb_2IOI.fasta", texto2 = "Cellular tumor antigen p53",
     nombre = "p53 humano vs p53 ratón")
]

for comp in comparaciones
    println("Procesando comparación: ", comp.nombre)
    secs1 = cargar_fasta_completo(comp.file1)
    secs2 = cargar_fasta_completo(comp.file2)

    header1, seq1 = buscar_secuencia_por_texto(secs1, comp.texto1)
    header2, seq2 = buscar_secuencia_por_texto(secs2, comp.texto2)

    

    alinear_y_mostrar(seq1, seq2, comp.nombre)

    dotplot(String(seq1), String(seq2), extraer_nombre(header1), extraer_nombre(header2), window=3)
end
