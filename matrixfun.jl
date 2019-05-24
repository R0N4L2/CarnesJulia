using Combinatorics
using DataValues
using LinearAlgebra
#Funcion generadora de matriz de equilibrio de peso por materia prima
function matriz_generador(M::Array{<:Any,2},carne::String,peso::Float64)
    P=M[2:end,:]
    replace!(x -> isna(x) ? 0 : x, P)
    n=size(P,1)
    R=P
    replace!(x -> x<0 ? 0 : x, R)
    Q=R
    a=findfirst(P[1,:].==0)
    for k=2:(n-3*(carne=="pA"))
        A=hcat(collect(combinations(2:n,k))...)'
        if carne=="pA"
            #restriccion que no se pueden combinar tomahauk con rib steak o
            #club steak
            C=A[(all!(trues(size(A,1)),A.!=9) .|all!(trues(size(A,1)),A.<10 .|A.>11)).&(all!(trues(size(A,1)),A.!=7) .|all!(trues(size(A,1)),A.<8 .| A.>11)),:]
        else
            C=A
        end
        if !isempty(C)
            i1=[unique(M[1,findall(any!(trues(size(P,2)-a+1),P[C[j,:],a:end]'.!=0)).+(a-1)])' for j=1:size(C,1)]
            i2=[sum(sum(Q[C[j,:],a:end],dims=2)) for j=1:size(C,1)]
            i3=[sum(Q[1,findall(any!(trues(a-1),in.(setdiff(M[1,findall(any!(trues(size(P,2)-a+1),P[C[j,:],a:end]'.!=0)).+(a-1)],[15273,15274])',M[1,1:a-1])))]) for j=1:size(C,1)]
            ie=findall(x->abs(x)<1e-5,i3-i2)
            if !isempty(ie)
                C=C[ie,:]
                i1=i1[ie]
                Z=vcat([hcat(Q[1,1:a-1]',sum(Q[C[j,:],a:end],dims=1)) for j=1:size(C,1)]...)
                ue=[findall(any!(trues(a-1),in.(i1[j],M[1,1:a-1]))) for j=1:size(i1,1)]
                row=findall(x->!isempty(x),ue)
                if !isempty(row)
                    for p=row for q=ue[p] Z[p,q]=0 end end
                    replace!(x -> x<0 ? 0 : x, Z)
                    R=vcat(R,Z)
                    R=unique(R[findall(x->abs(x)<1e-4,vcat(sum(R,dims=2).-peso...)),:],dims=1)
                end
            end
        end
    end
    f=findall(any!(trues(size(P,2)),P'.>0))
    R=R[:,f]
    b=P[1,f].==0
    n,m=size(R)
    D=zeros(n,m*n)
    for k=1:n,j=1:m
        D[k,j+(k-1)*m]=R[k,j]
    end
    D=hcat(D,diagm(0=>-peso*ones(n)))
    #A=vcat(trues(m*n),falses(n))'
    A=vcat(repeat(b,outer=n),falses(n))
    K=vcat(repeat(M[1,f],outer=n),-(1:n).-((carne=="pA")+2*(carne=="pB")+3*(carne=="bA"))
    *10^(2+ceil(log10(n))+(floor(log10(n))==ceil(log10(n)))))'
    ind=any!(trues(size(D,2)),D'.!=0)
    D=D[:,ind]
    A=A[ind]
    K=K[ind]
    B=vcat(M[1,f]',R)
    return K,D,A,B
end
