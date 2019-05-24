using PyCall
using ExcelReaders
using Combinatorics
using LinearAlgebra
using DataValues
using Statistics
using DataFrames
using XLSX
using Dates

esquema="INTEGRA"
option2="Igualdad en brazos y piernas"
#lam: peso que se le da a la restirccion de compra por cadenas
lam=.5
option3="Prioridad de demanda"
#pss: peso por sobrestock
pss=10
pss=float(pss)
#parametros de optimizacion
tiempo=360 #segundos
tiempo=float(tiempo)
gap=.1 #taza de aproximacion al optimo
tareas=10 #trabajos totales en paralelo
#crea reporte en xls
data_directory=pwd()
include(joinpath(data_directory,"matrixfun2.jl"))
MPcD=EjecutarQuery(join(["select DISTINCT replace(replace(a.DSC1,'DE RES CALIDAD (',''),')','') AS TIPO,a.ITM ",
"from ",esquema,".itt_i_maestros a, ",esquema,".itt_i_recetas b where upper(dsc1) like '%RES%' and upper(dl01) like '%MATERIA%'",
" and a.itm=b.itm and b.rscp>0"]))
MateriaPrimas=vec(convert(Array{Int64},MPcD.ITM))
MMP=hcat([reshape([j,GenRecetas(j,esquema)...],(4,1)) for j in MateriaPrimas]...)
MMP=MMP[:,findall(vec(prod(MMP[2:3,:],dims=1).>0))]
peso=EjecutarQuery(join(["select distinct itm, qnty/10000 as peso from ",esquema,".itt_i_recetas where itm in (",
string(convert(Array{Int64},MMP[1,:]))[2:end-1],")"]))#readxlsheet(joinpath(data_directory, "optimizacion_carnes3.xls"), "pesaje")
d2=sum([MMP[2,j][1]-1 for j=1:size(MMP,2)])+1
d1=sum([MMP[3,j][1] for j=1:size(MMP,2)])
pB=zeros(d1,d2)
c2=cumsum([MMP[2,j][1]-1 for j=1:size(MMP,2)]).+1
c1=cumsum([MMP[3,j][1] for j=1:size(MMP,2)])
for j=1:size(MMP,2)
    if j==1
        pB[1:c1[j],1]=MMP[4,j][1,:]
        pB[1:c1[j],2:c2[j]]=peso.PESO[MMP[1,j].==peso.ITM][1].*MMP[4,j][2:end,:]'
    else
        pB[(c1[j-1]+1):c1[j],1]=MMP[4,j][1,:]
        pB[(c1[j-1]+1):c1[j],(c2[j-1]+1):c2[j]]=peso.PESO[MMP[1,j].==peso.ITM][1].*MMP[4,j][2:end,:]'
    end
end
a=unique(pB[:,1])
pB=vcat([hcat(j,sum(pB[pB[:,1].==j,2:end],dims=1)) for j=a]...)
pB=vcat(hcat(-ones(2,1),hcat(vcat([0:(MMP[2,j]-2) for j=1:size(MMP,2)]...),vcat([MMP[1,j]*ones(MMP[2,j]-1) for j=1:size(MMP,2)]...))'),pB) #transformacion corte\[receta,materiaprima]
SE=convert(Array{Int},a)
dem=EjecutarQuery(join(["select a.drqj as fechapred,b.itm,sum(a.uorg*b.qnty*b.nvs)/100000000000 as dem_min,sum(a.uorg*b.qnty*(b.nof+10000))/100000000000",
" as sobre_oferta from ",esquema,".itt_i_forecast a , ",esquema,".itt_i_algoritmo b where b.kit=a.itm and b.itm in  (",string(SE)[2:end-1],
") group by b.itm,a.drqj order by a.drqj"]))
dem=convert(Matrix{Float64},dem)
MPcD=hcat(vec(convert(Array{String},MPcD.TIPO)),MateriaPrimas)
a=unique(dem[:,1])
OCC=hcat([reshape([j,OptimoCorte(dem[j.==dem[:,1],2:end],pB,option2,lam,MPcD,option3,tiempo,gap,tareas,pss,esquema)...],(3,1)) for j in a]...)
for n in a
    InserTable(convert(Array{Float64},hcat(OCC[2,OCC[1,:].==n][1][2:end,[4,2]],OCC[2,OCC[1,:].==n][1][2:end,5:end])),join([esquema,".ITT_O_WO"]),n)
end
OPT=vcat([OptimoProceso(OCC[3,OCC[1,:].==j][1][k,:],tiempo,gap,tareas,j,esquema) for j in unique(OCC[1,:]) for k=1:size(OCC[3,OCC[1,:].==j][1],1)]...)
if any(OPT[:,end].>=0)
    OPT=OPT[OPT[:,end].>=0,:]
    a=unique(OPT[:,1])
    OPT=convert(Array{Float64},OPT[:,[1,6,10]])
    for n in a
        InserTable(OPT[OPT[:,1].==n,2:3],join([esquema,".ITT_O_MPS"]),n)
    end
else
    print("ERROR:Problema Infactible")
end
