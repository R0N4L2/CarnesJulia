using Combinatorics
using DataValues
using LinearAlgebra
using ODBC
using Dates
using Cbc
using Ipopt
using Clp
using Juniper
using JuMP
#Funcion generadora de matriz de equilibrio de peso por materia prima
function GenMatriz(MateriaPrima::Int64,esquema::String)
    M=convert(Matrix{Float64},EjecutarQuery(join(["select qnrc/100,kitc,rscp/10000",
    " from ",esquema,".itt_i_recetas  where rscp>0 and itm=",string(MateriaPrima)])))
    if size(M,1)>0
        a=sort(unique(M[:,2]))
        b=sort(unique(M[:,1]))
        Z=vcat(a',zeros(size(b,1),size(a,1)))
        for (i,m) in enumerate(b) for (j,n) in enumerate(a) if any((M[:,1].==m).&(M[:,2].==n)) Z[i+1,j]=M[(M[:,1].==m).&(M[:,2].==n),3][1] end end end
        if size(Z,1)>2
            E=vcat([(Z[2,:]-Z[k,:])' for k=3:size(Z,1)]...)
            D=vcat(Z[1,:]',zeros(1,size(Z,2)),E)
            G=hcat([vcat(D[1,j],D[2:end,j].*(D[2:end,j].<0)) for j=1:size(D,2)]...)
            Z[2:end,:]+=G[2:end,:]
            Z=hcat(Z,vcat(G[1,:]',-G[2:end,:]))
            Z=Z[:,findall(sum(Z[2:end,:]',dims=2).>0)]
        end
    else
        Z=M
    end
    return Z
end
function GenRecetas(MateriaPrima::Int64,esquema::String)
    M=GenMatriz(MateriaPrima,esquema)
    molida=vec(convert(Matrix{Float64},EjecutarQuery(join(["select DISTINCT ITM from ",esquema,".itt_i_maestros where upper(dsc1)",
    " like '%MOLIDA%' and upper(dl01) like '%DESPIECE%' AND ITM IN ( SELECT KITC FROM ",esquema,".itt_i_recetas WHERE ITM=",string(MateriaPrima),")"]))))
    if size(M,1)>2
        z=findfirst(M[2,:].==0)
        if any(all!(trues(size(M,2)-z+1),M[1,z:end].!=molida'))
            rest=convert(Matrix{Int64},EjecutarQuery(join(["select qnrc/100+1,kitc from ",esquema,".itt_i_recetas where rscp<0 and itm=",string(MateriaPrima)])))
            R=M[2:end,:]
            n=size(R,1)
            Q=R
            p=minimum(sum(R,dims=2))
            P=maximum(sum(R,dims=2))
            #restriccion que no se pueden combinar tomahauk con rib steak o club steak
            if size(rest,1)>0
                a=unique(rest[:,1])
                b=unique(rest[:,2])
                K=zeros(size(a,1),size(b,1))
                for (i,m) in enumerate(a) for (j,n) in enumerate(b) K[i,j]=sum((rest[:,1].==m) .& (rest[:,2].==n)) end end
                K=convert(Array{Int64},hcat(a,sum(K,dims=2)))
                if any(K[:,2].>1)
                    s=K[K[:,2].>1,1][1]
                else
                    s=0
                end
                q=K[K[:,2].==1,1]
                t=size(K,1)
            else
                t=0
            end
            for k=2:(n-t) #negativos receta
                A=hcat(collect(combinations(2:n,k))...)'
                if t>0
                    #restriccion que no se pueden combinar tomahauk con rib steak o
                    #club steak, el bife de chorizo y todas las otras combinaciones del lomo de falda
                    C=A[all!(trues(size(A,1)),A.!=s).|all!(trues(size(A,1)),hcat([A.!=j for j=q]...)),:]
                else
                    C=A
                end
                if !isempty(C)
                    i1=[unique(M[1,findall(any!(trues(size(Q,2)-z+1),Q[C[j,:],z:end]'.>0)).+(z-1)])' for j=1:size(C,1)]
                    i2=[sum(Q[C[j,:],z:end]) for j=1:size(C,1)]
                    i3=[sum(Q[1,findall(any!(trues(z-1),in.(setdiff(M[1,findall(any!(trues(size(Q,2)-z+1),Q[C[j,:],z:end]'.>0)).+(z-1)],molida)',M[1,1:z-1])))]) for j=1:size(C,1)]
                    ie=findall(x->abs(x)<(P-p),i3-i2)
                    if !isempty(ie)
                        C=C[ie,:]
                        i1=i1[ie]
                        Z=vcat([hcat(Q[1,1:z-1]',sum(Q[C[j,:],z:end],dims=1)) for j=1:size(C,1)]...)
                        ue=[findall(any!(trues(z-1),in.(i1[j],M[1,1:z-1]))) for j=1:size(i1,1)]
                        row=findall(x->!isempty(x),ue)
                        if !isempty(row)
                            for p=row for q=ue[p] Z[p,q]=0 end end
                            if any(Z.>0) R=vcat(R,Z[any!(trues(size(Z,1)),Z.>0),:]) end
                            R=unique(R[vec((sum(R,dims=2).<=P).&(sum(R,dims=2).>=p)),:],dims=1)
                        end
                    end
                end
            end
            R=vcat(M[1,:]',R)
        else
            R=M
        end
    else
        R=M
    end
    return size(R)...,R
end
function EjecutarQuery(consulta::String)
    db=ODBC.DSN("JDE_ODBC_ORACLE","CONS_CARNICOS","C4rn1c0s")
    qry=ODBC.query(db,consulta)
    ODBC.disconnect!(db)
    return qry
end
function InserTable(X::Array{Float64},table::String,fechainsert::Float64)
    upmj=string(Date2JDE(today()))
    user="'OPTCARJL'"
    columnas="UPMJ,TDAY,USERM,DRQJ"
    fechainsert=strInt(fechainsert)
    if split(table,".")[2]=="ITT_O_MPS"
        columnas=join([columnas,"ITM,UORG"],",")
    else
        columnas=join(["UKID",columnas,"KIT,UORG"],",")
        contador=convert(Matrix{Any},EjecutarQuery(join(["select max(UKID) from ",table])))[1]
        if isequal(contador,missing)
            contador=0
        end
    end
    db=ODBC.DSN("JDE_ODBC_ORACLE","CONS_CARNICOS","C4rn1c0s")
    if split(table,".")[2]=="ITT_O_MPS"
        stmt=join(["INSERT INTO ",table," (",columnas,") VALUES ("])
        for row=1:size(X,1)
            tday=string(tiempoEntero(now()))
            stmt1=join([stmt,join(hcat(upmj,tday,user,fechainsert,strInt(X[row,1]),strFlt(X[row,2])),","),")"])
            try
                ODBC.execute!(db,stmt1)
            catch
                print(stmt1)
            end
        end
    else
        for row=2:size(X, 1)
            a=X[row,3:end].>0
            b=findall(a).+2
            columnas1=join([columnas,join([join(["ITM",j,",","UORG",j]) for j=1:1:sum(a)],",")],",")
            tday=string(tiempoEntero(now()))
            stmt=join(["INSERT INTO ",table," (", columnas1,") VALUES ("])
            stmt=join([stmt,join(hcat(strInt(contador+row-1),upmj,tday,user,fechainsert,strInt(X[row,1]),strInt(X[row,2]),hcat([hcat(strInt(X[1,j]),strFlt(X[row,j])) for j in b]...)),","),")"])
            try
                ODBC.execute!(db,stmt)
            catch
                print(stmt)
            end
        end
    end
    ODBC.disconnect!(db)
end
function strInt(n::Union{Int8,Int16,Int32,Int64,Int128,UInt32,UInt64,UInt8,UInt128,UInt16,Float16,Float32,Float64})
    return string(Int(round(n)))
end
function strFlt(n::Union{Float16,Float32,Float64})
    return string(Float64(n))
end
function Date2JDE(dat::Date)
    y=year(dat)
    return Int(round(((floor(y/100)-19)*100+mod(y,100))*1000+Dates.value(dat-Date(y-1,12, 31))))
end
function DtJDE2Date(dat::Int64)
    days=mod(dat,1000)
    y=mod(floor(dat/1000),100)+(floor(dat/100000)+19)*100
    dt=Dates.monthday(days)
    return replace(string(convert(Array{Int64},[dt[2:-1:1]...,y]))[2:end-1],", "=>"/")
end
function tiempoEntero(dt::DateTime)
    return Int((hour(dt)*10^2+minute(dt))*10^2+second(dt))
end
function OptimoCorte(dem::Array{Float64},pB::Array{Float64},
    option2::String,lam::Float64,MPcD::Array{Any},option3::String,
    tiempo::Float64,gap::Float64,tareas::Int64,pss::Float64,esquema::String)
    if size(dem,1)>0
        if size(dem,1)<size(pB,1)-2
            dems=hcat(pB[3:end,1],zeros(size(pB,1)-2,2))
            for j=1:size(dem,1) dems[dems[:,1].==dem[j,1],2:3]=dem[j,2:3] end
            dems=vcat(-ones(2,3),dems)
        end
        if option2=="Igualdad en brazos y piernas"
            a=findall(any!(trues(size(MPcD,1)),unique(pB[2,2:end])'.==MPcD[:,2]))
            v=hcat(MPcD[a,:],zeros(size(a,1),1))
            if any(v[:,1].=="PIERNA A") v[findfirst(v[:,1].=="PIERNA A"),3]=11 end
            if any(v[:,1].=="BRAZO A") v[findfirst(v[:,1].=="BRAZO A"),3]=12 end
            if any(v[:,1].=="PIERNA B") v[findfirst(v[:,1].=="PIERNA B"),3]=21 end
            if any(v[:,1].=="BRAZO B") v[findfirst(v[:,1].=="BRAZO B"),3]=22 end
            if any(v[:,1].=="PIERNA U") v[findfirst(v[:,1].=="PIERNA U"),3]=31 end
            if any(v[:,1].=="BRAZO U") v[findfirst(v[:,1].=="BRAZO U"),3]=32 end
            for i=1:3
                if (sum((v[:,3].>i*10).&(v[:,3].<(i+1)*10))!=2) & (any((v[:,3].>i*10).&(v[:,3].<(i+1)*10)))
                    v[(v[:,3].>i*10).&(v[:,3].<(i+1)*10),3]=0
                end
            end
            v=v[v[:,3].>0,:]
            M=zeros(Int(size(v,1)/2),size(pB,2)-1)
            for i=1:3
                if any((v[:,3].>i*10).&(v[:,3].<(i+1)*10))
                    M[Int(sum(v[:,3].<=i*10+2)/2),findall(v[v[:,3].==i*10+1,2].==pB[2,:]).-1]=ones(1,sum(v[v[:,3].==i*10+1,2].==pB[2,:]))
                    M[Int(sum(v[:,3].<=i*10+2)/2),findall(v[v[:,3].==i*10+2,2].==pB[2,:]).-1]=-ones(1,sum(v[v[:,3].==i*10+2,2].==pB[2,:]))
                end
            end
            dm=vcat(dems[3:end,2:3],zeros(size(M,1),2))
            Lam=vcat(ones(size(dems,1)-2,1),lam.*ones(size(M,1),1))
            M=vcat(pB[3:end,2:end],M)
        else
            dm=dems[3:end,2:3]
            M=pB[3:end,2:end]
            Lam=ones(size(dems,1)-2,1)
        end
        modele = Model(solver=CbcSolver(seconds=tiempo,ratioGap=gap,threads=tareas))
        @variable(modele,0<=x[1:(size(pB,2)-1)]<=10000,category=:Int)
        @variable(modele,0<=e<=10000)
        if option3=="Prioridad de demanda"
            fo=sum(mean(dems[3:end,2:3],dims=2).*pB[3:end,2:end],dims=1)./sum(pB[3:end,2:end],dims=1)
            @objective(modele, Min,pss*e-sum(fo*x))
        else
            @objective(modele, Min,e)
        end
        @constraint(modele,M*x+Lam.*e.>=dm[:,1])
        @constraint(modele,M*x-Lam.*e.<=dm[:,2])
        statuse=solve(modele)
        if statuse==:Optimal
            x_value=getvalue(x)
            e_value=getvalue(e)
            if any(x_value.>0)
                X=hcat(vcat(-1,x_value),pB')
                X=X[findall(x->abs(x)>1e-4,X[:,1]),:]
                X=X[:,vcat(1:3,findall(x->abs(x)>1e-4,vec(sum(X[2:end,4:end],dims=1))).+3)]
                MPcD=convert(Matrix{Any},EjecutarQuery(join(["select DISTINCT DSC1,ITM from ",esquema,".itt_i_maestros",
                " where itm in (",string(convert(Array{Int64},X[2:end,3]))[2:end-1],")"])))
                X=hcat(vcat(-1,MPcD[[findfirst(MPcD[:,2].==j) for j in X[2:end,3]],1]),X)
                SEcD=convert(Matrix{Any},EjecutarQuery(join(["select DISTINCT DSC1,ITM from ",esquema,".itt_i_maestros",
                " where itm in (",string(convert(Array{Int64},X[1,5:end]))[2:end-1],")"])))
                X=vcat(hcat(-ones(1,4),reshape(SEcD[[findfirst(SEcD[:,2].==j) for j in X[1,5:end]],1],1,size(X,2)-4)),X)
                Y=vcat([hcat(reshape(X[1:2,j],1,2),sum(X[3:end,2].*X[3:end,j])) for j=5:size(X,2)]...)
            else
                X="No tiene cortes optimos"
                Y=X
            end
        else
            X="Problema Infactible"
            Y=X
        end
    else
        X="no cortar"
        Y=X
    end
    return X,Y
end
function OptimoProceso(Y::Array{Any},tiempo::Float64,gap::Float64,tareas::Int64,fecha::Float64,esquema::String)
    dem=EjecutarQuery(join(["select c.dsc1, a.itm,a.uorg*b.nvs/(b.qnty*1000) as dem_min,a.uorg*(b.nof+10000)/(b.qnty*1000) ",
    "as sobre_oferta,10000/b.qnty as proporcion from ",esquema,".itt_i_forecast a,",
    esquema,".itt_i_algoritmo b,",esquema,".itt_i_maestros c where b.kit=a.itm and b.itm=",
    strInt(Y[2])," and a.drqj=",strInt(fecha)," and a.itm=c.itm"]))
    if size(dem,1)>0
        dem=convert(Matrix{Any},dem)
        dm=convert(Matrix{Float64},dem[:,3:5])
        modele=Model(solver=CbcSolver(seconds=tiempo,ratioGap=gap,threads=tareas))
        @variable(modele,0<=x[1:size(dm,1)]<=10000,category=:Int)
        @variable(modele,0<=e<=10000)
        M=Matrix(vcat(Diagonal(ones(size(dm,1))),reshape(dm[:,3],1,size(dm,1))))
        @objective(modele,Min,e)
        @constraint(modele,M*x.+e.>=vcat(dm[:,1],Y[3]))
        @constraint(modele,M*x.-e.<=vcat(dm[:,2],Y[3]))
        statuse=solve(modele)
        if statuse==:Optimal
            x_value=getvalue(x)
            e_value=getvalue(e)
            dem=hcat(repeat(hcat(fecha,reshape(Y,1,prod(size(Y)))),size(dem,1),1),dem,x_value)
        else
            dem=hcat(repeat(hcat(fecha,reshape(Y,1,prod(size(Y)))),size(dem,1),1),dem,-ones(size(dem,1)))
        end
    else
        dem=hcat(fecha,reshape(Y,1,prod(size(Y))),-ones(1,6))
    end
    return dem
end
