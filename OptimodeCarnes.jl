using PyCall
using ExcelReaders
using Combinatorics
using LinearAlgebra
using DataValues
using Statistics
import DataFrames, XLSX
#data_directory = joinpath(@__DIR__, "..", "data")
data_directory ="C:\\Users\\ronal\\OneDrive\\Documentos\\carnes\\"
include(joinpath(data_directory,"matrixfun.jl"))
pA=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "Receta piernaA")
pB=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "Receta piernaB")
bA=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "Receta brazoA")
peso=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "pesaje")
demstr=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "Demanda")
fostr=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "NivelServicio")
dic=readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "Diccionario")
files=readdir(data_directory)
pbA=peso[2,2]#peso del brazo tipo A
ppA=peso[3,2]#peso de la pierna tipo A
ppB=peso[4,2]#peso de la pierna tipo B
dem=vcat([hcat(j,sum(demstr[demstr[:,4].==j,3])) for j=unique(demstr[2:end,4])]...)
q=.!any!(trues(size(dic,1)-2),in.(dic[2:end-1,1],dem[:,1]'))
if any(q)
    dem=vcat(dem,hcat(dic[findall(q).+1],zeros(sum(q))))
end
a=unique(fostr[2:end,7])
fo=vcat([hcat(j,sum(fostr[fostr[:,7].==j,6])/sum(mean(fostr[fostr[:,7].==j,3:4],dims=2).+.5)) for j=a]...)
fp=vcat([hcat(j,mean(fostr[fostr[:,7].==j,3])) for j=a]...)
fq=vcat([hcat(j,mean(fostr[fostr[:,7].==j,4].+1)) for j=a]...)
ti1,MpA,rel1,C1=matriz_generador(pA,"pA",ppA)
ti2,MpB,rel2,C2=matriz_generador(pB,"pB",ppB)
ti3,MbA,rel3,C3=matriz_generador(bA,"bA",pbA)
s0,d0=size(MpA)
s1,d1=size(MpB)
s2,d2=size(MbA)
Meq=vcat(hcat(MpA,zeros(s0,d1+d2)),hcat(zeros(s1,d0),MpB,zeros(s1,d2)),hcat(zeros(s2,d0+d1),MbA)) #restircciones de equilibrio de pesos
Mpb=hcat(zeros(1,d0-s0),ones(1,s0),zeros(1,d1),zeros(1,d2-s2),-ones(1,s2)) #Matriz de diferencias entre priernas y brazos, para compra por cadenas
tpA=findall(ti1.<0)
tpB=findall(ti2.<0).+d0
tbA=findall(ti3.<0).+(d0+d1)
tipo=vcat(ti1,ti2,ti3)
f1=zeros(length(tipo))
for j=1:length(tipo)
    if in(tipo[j],fo[:,1])
        f1[j]=sum(fo[in.(tipo[j],fo[:,1]),2])+1
    elseif in(j,tpA)
        f1[j]=-ppA-10
    elseif in(j,tpB)
        f1[j]=-ppB-5
    elseif in(j,tbA)
        f1[j]=-pbA-2
    else
        f1[j]=.5
    end
end
if any(dem[:,2].>0)
    demp=dem[dem[:,2].>0,:]
    Mi=zeros(size(demp,1),length(tipo))
    bi=demp[:,2]
    fi=zeros(size(bi,1))
    fj=zeros(size(bi,1))
    for j=1:size(demp,1)
        fi[j]=sum(fp[fp[:,1].==demp[j,1],2])
        fj[j]=sum(fq[fq[:,1].==demp[j,1],2])
        a=findall(in.(demp[j,1],tipo))
        Z=sum(unique(Meq[:,a],dims=1),dims=1)
        for k=1:length(a)
            Mi[j,a[k]]=Z[k]
        end
    end
end
if any(dem[:,2].==0)
    dem0=dem[dem[:,2].==0,:]
    Meq0=zeros(size(dem0,1),length(tipo))
    beq0=dem0[:,2]
    for j=1:size(dem0,1)
        a=findall(in.(dem0[j,1],tipo))
        Z=sum(unique(Meq[:,a],dims=1),dims=1)
        for k=1:length(a)
            Meq0[j,a[k]]=Z[k]
        end
    end
end
include(joinpath(data_directory,"OLP.jl"))
restri="igualdad"
Rel=vcat(rel1,rel2,rel3)
status,x_value,e=LP(f1,Mi,bi,fi,fj,Meq0,Meq,Mpb,.5,restri)
if status!=:Optimal
    restri="desigualdad"
end
status,x_value,e=LP(f1,Mi,bi,fi,fj,Meq0,Meq,Mpb,.5,restri)
# elimina variables y restricciones innecesarias para el problema lineal,
# al simplificar el problema lineal con menos variables, el problema entero
# tambien tendra menos variables, por lo que el arbol de ramificacion y corte
# será mucho mas pequeo y manejable, por lo que se demorara mucho menos en
# encontrar la solucion entera factible optima
if status==:Optimal
        stop=true
        while stop
            if status==:Optimal
                idx=(x_value.<10000).&((x_value.>0).|(tipo.>0))
                if !all(idx)&(length(idx)>2)
                        di=any!(trues(sum(tipo.<0)),in.(tipo[tipo.<0],tipo[(tipo.<0) .&idx]'))
                        Mi=Mi[:,idx]
                        Meq=Meq[di,idx]
                        Meq0=Meq0[:,idx]
                        f1=f1[idx]
                        tipo=tipo[idx]
                        Mpb=Mpb[:,idx]
                        idx=any!(trues(size(Meq,2)),Meq'.!=0)
                    end
                else
                    stop=false
                end
            else
                stop=false
            end
        end
    ## resultado lineal
    hc=hcat(tipo[x_value.>0],x_value[x_value.>0],vcat(map(x->hcat(x[1],x[2]),findall(Meq[:,x_value.>0].!=0))...))
    df1=DataFrames.DataFrame(CODIGO=hc[:,1],Valor=hc[:,2],fila=hc[:,3],columna=hc[:,4])
    Conf0pA=C1[vcat(true,any!(trues(sum(ti1.<0)),in.(ti1[ti1.<0],tipo[tipo.<0]'))),:]
    Conf0pA=hcat(Conf0pA,vcat(0,x_value[any!(trues(size(tipo,1)),in.(tipo,ti1[ti1.<0]'))]))
    Conf0pA=Conf0pA[vcat(1,findall(Conf0pA[2:end,end].>0 .*Conf0pA[2:end,end].<10000).+1),:]
    Conf0pA=Conf0pA[:,any!(trues(size(Conf0pA,2)),Conf0pA[2:end,:]'.>0)]
    a=sort(unique(Conf0pA[1,:]))
    Conf0pA=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0pA[2:end,Conf0pA[1,:].==j],dims=2)) for j=a]...)
    Conf0pA[1,1]="A-Cantidad de pierna A"
    Conf0pB=C2[vcat(true,any!(trues(sum(ti2.<0)),in.(ti2[ti2.<0],tipo[tipo.<0]'))),:]
    Conf0pB=hcat(Conf0pB,vcat(0,x_value[any!(trues(size(tipo,1)),in.(tipo,ti2[ti2.<0]'))]))
    Conf0pB=Conf0pB[vcat(1,findall(Conf0pB[2:end,end].>0 .*Conf0pB[2:end,end].<10000).+1),:]
    Conf0pB=Conf0pB[:,any!(trues(size(Conf0pB,2)),Conf0pB[2:end,:]'.>0)]
    a=sort(unique(Conf0pB[1,:]))
    Conf0pB=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0pB[2:end,Conf0pB[1,:].==j],dims=2)) for j=a]...)
    Conf0pB[1,1]="A-Cantidad de pierna B"
    Conf0bA=C3[vcat(true,any!(trues(sum(ti3.<0)),in.(ti3[ti3.<0],tipo[tipo.<0]'))),:]
    Conf0bA=hcat(Conf0bA,vcat(0,x_value[any!(trues(size(tipo,1)),in.(tipo,ti3[ti3.<0]'))]))
    Conf0bA=Conf0bA[vcat(1,findall(Conf0bA[2:end,end].>0 .*Conf0bA[2:end,end].<10000).+1),:]
    Conf0bA=Conf0bA[:,any!(trues(size(Conf0bA,2)),Conf0bA[2:end,:]'.>0)]
    a=sort(unique(Conf0bA[1,:]))
    Conf0bA=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0bA[2:end,Conf0bA[1,:].==j],dims=2)) for j=a]...)
    Conf0bA[1,1]="A-Cantidad de brazo A"
    df2=DataFrames.DataFrame(Nombre_producto=vcat([dic[j.==dic[:,1],2] for j=vcat(demp[:,1],dem0[:,1])]),CODIGO=vcat(demp[:,1],dem0[:,1]),
    Resultado=vcat(Mi*x_value,Meq0*x_value),demanda=vcat(bi,zeros(size(Meq0,1))),
    cota_inf=vcat(fi,zeros(size(Meq0,1))),cota_sup=vcat(fj,zeros(size(Meq0,1))))
    df3=DataFrames.DataFrame(Dict(Conf0pA[1,j]=>Conf0pA[2:end,j] for j=1:size(Conf0pA,2)))
    df4=DataFrames.DataFrame(Dict(Conf0pB[1,j]=>Conf0pB[2:end,j] for j=1:size(Conf0pB,2)))
    df5=DataFrames.DataFrame(Dict(Conf0bA[1,j]=>Conf0bA[2:end,j] for j=1:size(Conf0bA,2)))
    if in("resultado_lineal.xlsx",files)
        rm(joinpath(data_directory,"resultado_lineal.xlsx"),force=true)
    end
    XLSX.writetable(joinpath(data_directory,"resultado_lineal.xlsx"),Valor_codigo=(DataFrames.columns(df1), DataFrames.names(df1)),
    demanda=(DataFrames.columns(df2),DataFrames.names(df2)),cortes_piernaA=(DataFrames.columns(df3),DataFrames.names(df3)),
    cortes_piernaB=(DataFrames.columns(df4),DataFrames.names(df4)),cortes_brazoA=(DataFrames.columns(df5),DataFrames.names(df5)))
    print("Optimizacion lineal realizada")
    ## problema entero pequenho
    if false
        statuse,xe,e=IP(Rel,f1,Mi,bi,fi,Meq0,Meq,Mpb,.5,restri)
        if statuse==:Optimal
            hc=hcat(tipo[xe.>0],xe[xe.>0],vcat(map(x->hcat(x[1],x[2]),findall(Meq[:,xe.>0].!=0))...))
            df1=DataFrames.DataFrame(CODIGO=hc[:,1],Valor=hc[:,2],fila=hc[:,3],columna=hc[:,4])
            hc=hc[Rel[xe.>0]]
            AA=repeat(Rel.|(tipo.<0),outer=size(Meq,1)).*(Meq.>0-Meq.<0).*repeat(xe',outer=size(Meq,1))
            yr=sum(repeat(Rel,outer=size(Meq,1)).*AA./repeat(sum(AA.*(AA.<0),dims=2),outer=(1,size(Meq,2))),dims=1)
            Conf0pA=C1[vcat(true,any!(trues(sum(ti1.<0)),in.(ti1[ti1.<0],tipo[tipo.<0]'))),:]
            Conf0pA=hcat(Conf0pA,vcat(0,xe[any!(trues(size(tipo,1)),in.(tipo,ti1[ti1.<0]'))]))
            Conf0pA=Conf0pA[vcat(1,findall(Conf0pA[2:end,end].>0 .*Conf0pA[2:end,end].<10000).+1),:]
            Conf0pA=Conf0pA[:,any!(trues(size(Conf0pA,2)),Conf0pA[2:end,:]'.>0)]
            a=sort(unique(Conf0pA[1,:]))
            Conf0pA=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0pA[2:end,Conf0pA[1,:].==j],dims=2)) for j=a]...)
            Conf0pA[1,1]="A-Cantidad de pierna A"
            Conf0pB=C2[vcat(true,any!(trues(sum(ti2.<0)),in.(ti2[ti2.<0],tipo[tipo.<0]'))),:]
            Conf0pB=hcat(Conf0pB,vcat(0,xe[any!(trues(size(tipo,1)),in.(tipo,ti2[ti2.<0]'))]))
            Conf0pB=Conf0pB[vcat(1,findall(Conf0pB[2:end,end].>0 .*Conf0pB[2:end,end].<10000).+1),:]
            Conf0pB=Conf0pB[:,any!(trues(size(Conf0pB,2)),Conf0pB[2:end,:]'.>0)]
            a=sort(unique(Conf0pB[1,:]))
            Conf0pB=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0pB[2:end,Conf0pB[1,:].==j],dims=2)) for j=a]...)
            Conf0pB[1,1]="A-Cantidad de pierna B"
            Conf0bA=C3[vcat(true,any!(trues(sum(ti3.<0)),in.(ti3[ti3.<0],tipo[tipo.<0]'))),:]
            Conf0bA=hcat(Conf0bA,vcat(0,xe[any!(trues(size(tipo,1)),in.(tipo,ti3[ti3.<0]'))]))
            Conf0bA=Conf0bA[vcat(1,findall(Conf0bA[2:end,end].>0 .*Conf0bA[2:end,end].<10000).+1),:]
            Conf0bA=Conf0bA[:,any!(trues(size(Conf0bA,2)),Conf0bA[2:end,:]'.>0)]
            a=sort(unique(Conf0bA[1,:]))
            Conf0bA=hcat([vcat(dic[j.==dic[:,1],2],sum(Conf0bA[2:end,Conf0bA[1,:].==j],dims=2)) for j=a]...)
            Conf0bA[1,1]="A-Cantidad de brazo A"
            hc=hcat(tipo[xe.>0],x_value[xe.>0],vcat(map(x->hcat(x[1],x[2]),findall(Meq[:,xe.>0].!=0))...))
            df1=DataFrames.DataFrame(CODIGO=hc[:,1],Valor=hc[:,2],fila=hc[:,3],columna=hc[:,4])
            hc=hc[Rel[xe.>0]]
            df2=DataFrames.DataFrame(Nombre_producto=vcat([dic[j.==dic[:,1],2] for j=vcat(demp[:,1],dem0[:,1])]),CODIGO=vcat(demp[:,1],dem0[:,1]),
            Resultado=vcat(Mi*xe,Meq0*xe),demanda=vcat(bi,zeros(size(Meq0,1))),
            conta=vcat(fi,zeros(size(Meq0,1))))
            df3=DataFrames.DataFrame(Dict(Conf0pA[1,j]=>Conf0pA[2:end,j] for j=1:size(Conf0pA,2)))
            df4=DataFrames.DataFrame(Dict(Conf0pB[1,j]=>Conf0pB[2:end,j] for j=1:size(Conf0pB,2)))
            df5=DataFrames.DataFrame(Dict(Conf0bA[1,j]=>Conf0bA[2:end,j] for j=1:size(Conf0bA,2)))
            if in("resultado_entero.xlsx",files)
                rm(joinpath(data_directory,"resultado_entero.xlsx"),force=true)
            end
            XLSX.writetable(joinpath(data_directory,"resultado_entero.xlsx"),Valor_codigo=(DataFrames.columns(df1), DataFrames.names(df1)),
            demanda=(DataFrames.columns(df2),DataFrames.names(df2)),cortes_piernaA=(DataFrames.columns(df3),DataFrames.names(df3)),
            cortes_piernaB=(DataFrames.columns(df4),DataFrames.names(df4)),cortes_brazoA=(DataFrames.columns(df5),DataFrames.names(df5)))
        else
            print("Error Optimo entero no encontrado")
        end
    end
else
    print("Error Optimo lineal no encontrado")
end
