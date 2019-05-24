using Pkg
Pkg.add("FFTW")
Pkg.add("Combinatorics")
Pkg.add("GSL")


isdir(Pkg.dir("Combinatorics"))
"Combinatorics" ∈ keys(Pkg.installed())
