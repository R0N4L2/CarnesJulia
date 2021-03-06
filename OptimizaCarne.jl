using Pkg
try
    if string(Pkg.installed()["PyCall"])!="1.91.1"
        Pkg.pin(PackageSpec(name="PyCall", version="1.91.1"))
    end
    using PyCall
catch
    Pkg.add("PyCall")
    Pkg.pin(PackageSpec(name="PyCall", version="1.91.1"))
    using PyCall
end
try
    if string(Pkg.installed()["Combinatorics"])!="0.7.0"
        Pkg.pin(PackageSpec(name="Combinatorics", version="0.7.0"))
    end
    using Combinatorics
catch
    Pkg.add("Combinatorics")
    Pkg.pin(PackageSpec(name="Combinatorics", version="0.7.0"))
    using Combinatorics
end
try
    using LinearAlgebra
catch
    Pkg.add("LinearAlgebra")
    using LinearAlgebra
end
try
    if string(Pkg.installed()["DataValues"])!="0.4.7"
        Pkg.pin(PackageSpec(name="DataValues",version="0.4.7"))
    end
    using DataValues
catch
    Pkg.add("DataValues")
    Pkg.pin(PackageSpec(name="DataValues",version="0.4.7"))
    using DataValues
end
try
    using Statistics
catch
    Pkg.add("Statistics")
    using Statistics
end
try
    if string(Pkg.installed()["DataFrames"])!="0.17.1"
        Pkg.pin(PackageSpec(name="DataFrames",version="0.17.1"))
    end
    using DataFrames
catch
    Pkg.add("DataFrames")
    Pkg.pin(PackageSpec(name="DataFrames",version="0.17.1"))
    using DataFrames
end
try
    using Dates
catch
    Pkg.add("Dates")
    using Dates
end
try
    if string(Pkg.installed()["ODBC"])!="0.8.1"
        Pkg.pin(PackageSpec(name="ODBC",version="0.8.1"))
    end
    using ODBC
catch
    Pkg.add("ODBC")
    Pkg.pin(PackageSpec(name="ODBC",version="0.8.1"))
    using ODBC
end
try
    if string(Pkg.installed()["Cbc"])!="0.6.0"
        Pkg.pin(PackageSpec(name="Cbc",version="0.6.0"))
    end
    using Cbc
catch
    Pkg.add("Cbc")
    Pkg.pin(PackageSpec(name="Cbc",version="0.6.0"))
    using Cbc
end
try
    if string(Pkg.installed()["JuMP"])!="0.18.5"
        Pkg.pin(PackageSpec(name="JuMP", version="0.18.5"))
    end
    using JuMP
catch
    using Pkg
    Pkg.add("JuMP")
    Pkg.pin(PackageSpec(name="JuMP", version="0.18.5"))
    using JuMP
end
function OptimoCarnes(esquema::String,option2::Bool,tiempo::Float64,gap::Float64,tareas::Int64,conn::Array{String})
    MPcD=EjecutarQuery(join(["select DISTINCT replace(replace(a.DSC1,'DE RES CALIDAD (',''),')','') AS TIPO,a.ITM ",
    "from ",esquema,".itt_i_maestros a, ",esquema,".itt_i_recetas b where upper(dsc1) like '%RES%' and upper(dl01) like '%MATERIA%'",
    " and a.itm=b.itm and b.rscp>0"]),conn)
    MateriaPrimas=vec(convert(Array{Int64},MPcD.ITM))
    MMP=hcat([reshape([j,GenRecetas(j,esquema,conn)...],(4,1)) for j in MateriaPrimas]...)
    MMP=MMP[:,findall(vec(prod(MMP[2:3,:],dims=1).>0))]
    peso=EjecutarQuery(join(["select distinct itm, qnty/10000 as peso from ",esquema,".itt_i_recetas where itm in (",
    string(convert(Array{Int64},MMP[1,:]))[2:end-1],")"]),conn)
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
    dem=EjecutarQuery(join(["select a.drqj as fechapred,b.itm,sum(a.uorg*b.qnty*b.nvs)*power(10,-12) as dem_min,sum(a.uorg*b.qnty*(b.nof+power(10,4)))*power(10,-12) as sobre_oferta from ",esquema,".itt_i_forecast a , ",esquema,".itt_i_algoritmo b where b.kit=a.itm and b.itm in  (",string(SE)[2:end-1],") and a.uorg>0 group by b.itm,a.drqj order by a.drqj"]),conn)
    if size(dem,1)>0
        dem=convert(Matrix{Float64},dem)
        MPcD=hcat(vec(convert(Array{String},MPcD.TIPO)),MateriaPrimas)
        a=unique(dem[:,1])
        OCC=hcat([reshape([j,OptimoCorte(dem[j.==dem[:,1],2:end],pB,option2,MPcD,tiempo,gap,tareas,esquema,conn)...],(3,1)) for j in a]...)
        for n in a
            if typeof(OCC[2,OCC[1,:].==n][1][1])<:Number
                InserTable(convert(Array{Float64},hcat(OCC[2,OCC[1,:].==n][1][2:end,[4,2]],OCC[2,OCC[1,:].==n][1][2:end,5:end])),join([esquema,".ITT_O_WO"]),n,conn)
            else
                println(join([DtJDE2Date(Int(n))," : ",OCC[2,OCC[1,:].==n][1],"\n"]))
            end
        end
        OCC=OCC[:,findall([typeof(OCC[2,OCC[1,:].==n][1][1])<:Number for n in a])]
        if !isempty(OCC)||size(OCC,1)>0
            OPT=vcat([OptimoProceso(OCC[3,OCC[1,:].==j][1][k,:],tiempo,gap,tareas,j,esquema,conn) for j in unique(OCC[1,:]) for k=1:size(OCC[3,OCC[1,:].==j][1],1)]...)
            if (!isempty(OPT)||size(OPT,1)>0) && any(isa.(OPT[:,end],Number))
                OPT1=OPT[isa.(OPT[:,end],Number),[1,6,end]]
                OPT1=convert(Array{Float64},OPT1)
                a=unique(OPT1[:,1])
                for n in a
                    InserTable(OPT1[OPT1[:,1].==n,2:3],join([esquema,".ITT_O_MPS"]),n,conn)
                end
            end
            if (!isempty(OPT)||size(OPT,1)>0) && any(isa.(OPT[:,end],String))
                OPT1=OPT[isa.(OPT[:,end],String),[1,2,3,end]]
                for n=1:size(OPT1,1)
                    println(join([DtJDE2Date(Int(OPT1[n,1]))," : ",OPT1[n,4]," de ",OPT1[n,2]," (codigo semielaborado:",strInt(OPT1[n,3]),")\n"]))
                end
            end
            if isempty(OPT)||size(OPT,1)==0
                println("ERROR:No existe factibilidad en ningun dia\n")
            end
        else
            println("ERROR:No existe factibilidad en ningun dia\n")
        end
    else
        println("ERROR:No existe demanda\n")
    end
end
function GenMatriz(MateriaPrima::Int64,esquema::String,conn::Array{String})
    M=convert(Matrix{Float64},EjecutarQuery(join(["select qnrc/100,kitc,rscp*power(10,-4) from ",esquema,".itt_i_recetas  where rscp>0 and itm=",string(MateriaPrima)]),conn))
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
function GenRecetas(MateriaPrima::Int64,esquema::String,conn::Array{String})
    M=GenMatriz(MateriaPrima,esquema,conn)
    molida=vec(convert(Matrix{Float64},EjecutarQuery(join(["select DISTINCT ITM from ",esquema,".itt_i_maestros where upper(dsc1)",
    " like '%MOLIDA%' and upper(dl01) like '%DESPIECE%' AND ITM IN ( SELECT KITC FROM ",esquema,".itt_i_recetas WHERE ITM=",string(MateriaPrima),")"]),conn)))
    if size(M,1)>2
        z=findfirst(M[2,:].==0)
        if any(all!(trues(size(M,2)-z+1),M[1,z:end].!=molida'))
            rest=convert(Matrix{Int64},EjecutarQuery(join(["select qnrc/100+1,kitc,rscp from ",esquema,".itt_i_recetas where rscp<0 and itm=",string(MateriaPrima)]),conn))
            R=M[2:end,:]
            n=size(R,1)
            Q=R
            p=minimum(sum(R,dims=2))
            P=maximum(sum(R,dims=2))
            #restriccion que no se pueden combinar tomahauk con rib steak o club steak
            if size(rest,1)>0
                if sum(minimum(rest[:,3])==rest[:,3])<2
                    t=size(rest,1)
                    s=rest[rest[:,3].==minimum(rest[:,3]),1][1]
                    q=rest[rest[:,3].!=minimum(rest[:,3]),1]
                else
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
                end
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
function EjecutarQuery(consulta::String,conn::Array{String})
    dsn,username,password=conn
    db=ODBC.DSN(dsn,username,password)
    qry=ODBC.query(db,consulta)
    ODBC.disconnect!(db)
    return qry
end
function InserTable(X::Array{Float64},table::String,fechainsert::Float64,conn::Array{String})
    upmj=string(Date2JDE(today()))
    user="'OPTCARJL'"
    columnas="UPMJ,TDAY,USERM,DRQJ"
    fechainsert=strInt(fechainsert)
    if split(table,".")[2]=="ITT_O_MPS"
        columnas=join([columnas,"ITM,UORG"],",")
    else
        columnas=join(["UKID",columnas,"KIT,UORG"],",")
        contador=convert(Matrix{Any},EjecutarQuery(join(["select max(UKID) from ",table]),conn))[1]
        if isequal(contador,missing)
            contador=0
        end
    end
    dsn,username,password=conn
    db=ODBC.DSN(dsn,username,password)
    tday=string(tiempoEntero(now()))
    if split(table,".")[2]=="ITT_O_MPS"
        stmt=join(["INSERT INTO ",table," (",columnas,") VALUES ("])
        for row=1:size(X,1)
            stmt1=join([stmt,join(hcat(upmj,tday,user,fechainsert,strInt(X[row,1]),strFlt(X[row,2]*10^4)),","),")"])
            try
                ODBC.execute!(db,stmt1)
            catch
                println(join([repeat("#",max(43,length(stmt1))),"\n"]))
                println("No puede insertar la siguiente linea:\n")
                println(join([stmt1,"\n"]))
                println("Revise si esa data ya esta en base de datos\n")
                println(join([repeat("-",max(43,length(stmt1))),"\n"]))
            end
        end
    else
        for row=2:size(X, 1)
            a=X[row,3:end].>0
            b=findall(a).+2
            columnas1=join([columnas,join([join(["ITM",j,",","UORG",j]) for j=1:1:sum(a)],",")],",")
            stmt=join(["INSERT INTO ",table," (", columnas1,") VALUES ("])
            stmt=join([stmt,join(hcat(strInt(contador+row-1),upmj,tday,user,fechainsert,strInt(X[row,1]),strInt(X[row,2]*10^4),hcat([hcat(strInt(X[1,j]),strFlt(X[row,j]*10^4)) for j in b]...)),","),")"])
            try
                ODBC.execute!(db,stmt)
            catch
                println(join([repeat("#",max(43,length(stmt))),"\n"]))
                println("No puede hacer el siguiente insert:\n")
                println(join([stmt,"\n"]))
                println("Revise si esa data ya esta en base de datos\n")
                println(join([repeat("-",max(43,length(stmt))),"\n"]))
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
    option2::Bool,MPcD::Array{Any},
    tiempo::Float64,gap::Float64,tareas::Int64,esquema::String,conn::Array{String})
    if size(dem,1)>0
        dems=hcat(pB[3:end,1],zeros(size(pB,1)-2,2))
        for j=1:size(dem,1)
            dems[dems[:,1].==dem[j,1],2:3]=dem[j,2:3]
        end
        if option2
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
            N=zeros(Int(size(v,1)/2),size(pB,2)-1)
            for i=1:3
                if any((v[:,3].>i*10).&(v[:,3].<(i+1)*10))
                    N[Int(sum(v[:,3].<=i*10+2)/2),findall(v[v[:,3].==i*10+1,2].==pB[2,:]).-1]=ones(1,sum(v[v[:,3].==i*10+1,2].==pB[2,:]))
                    N[Int(sum(v[:,3].<=i*10+2)/2),findall(v[v[:,3].==i*10+2,2].==pB[2,:]).-1]=-ones(1,sum(v[v[:,3].==i*10+2,2].==pB[2,:]))
                end
            end
        end
        dm=dems[:,2:3]
        M=pB[3:end,2:end]
        fo=(M'*mean(dm,dims=2))
        bl=vec((fo.>0).&(sum(M,dims=1)'.>0))
        M=M[:,bl]
        if option2
            N=N[:,bl]
        end
        fo=fo[bl]
        fo=sum(M,dims=1).*fo'.^-1
        modele= Model(solver=CbcSolver(seconds=tiempo,ratioGap=gap,threads=tareas))
        @variable(modele,0<=x[1:(size(M,2))]<=10^5,category=:Int)
        @variable(modele,0<=e<=10^5)
        if option2
            @variable(modele,0<=e1<=10^5)
        end
        @objective(modele,Min,sum(fo*x))
        @constraint(modele,M*x.>=dm[:,1])
        @constraint(modele,M*x.-e.<=dm[:,2])
        @constraint(modele,sum(M,dims=1)*x.>=.9*sum(dm[:,1:2])/2)
        @constraint(modele,sum(M,dims=1)*x.-e.<=1.08*sum(dm[:,1:2])/2)
        if option2
            @constraint(modele,N*x.-e1.<=0)
            @constraint(modele,N*x.+e1.>=0)
        end
        statuse=solve(modele)
        if statuse==:Optimal
            x_value=getvalue(x)
            e_value=getvalue(e)
            if option2
                e1_value=getvalue(e1)
            end
            if any(x_value.>0)
                X=hcat(vcat(-1,x_value),pB')
                X=X[findall(x->abs(x)>1e-4,X[:,1]),:]
                X=X[:,vcat(1:3,findall(x->abs(x)>1e-4,vec(sum(X[2:end,4:end],dims=1))).+3)]
                MPcD=convert(Matrix{Any},EjecutarQuery(join(["select DISTINCT DSC1,ITM from ",esquema,".itt_i_maestros",
                " where itm in (",string(convert(Array{Int64},X[2:end,3]))[2:end-1],")"]),conn))
                X=hcat(vcat(-1,MPcD[[findfirst(MPcD[:,2].==j) for j in X[2:end,3]],1]),X)
                SEcD=convert(Matrix{Any},EjecutarQuery(join(["select DISTINCT DSC1,ITM from ",esquema,".itt_i_maestros",
                " where itm in (",string(convert(Array{Int64},X[1,5:end]))[2:end-1],")"]),conn))
                X=vcat(hcat(-ones(1,4),reshape(SEcD[[findfirst(SEcD[:,2].==j) for j in X[1,5:end]],1],1,size(X,2)-4)),X)
                Y=vcat([hcat(reshape(X[1:2,j],1,2),sum(X[3:end,2].*X[3:end,j])) for j=5:size(X,2)]...)
            else
                X="No tiene cortes optimos"
                Y=X
                println(join([Y,"\n"]))
            end
        else
            X="Problema Infactible de corte"
            Y=X
            println(join([Y,"\n"]))
        end
    else
        X="no cortar"
        Y=X
        println(join([Y,"\n"]))
    end
    return X,Y
end
function OptimoProceso(Y::Array{Any},tiempo::Float64,gap::Float64,tareas::Int64,fecha::Float64,esquema::String,conn::Array{String})
    dem=EjecutarQuery(join(["select c.dsc1, a.itm,power(10,-4)*a.uorg*b.nvs/b.qnty as dem_min,power(10,-4)*a.uorg*(b.nof+power(10,4))/b.qnty as sobre_oferta,power(10,4)/b.qnty as proporcion from ",esquema,".itt_i_forecast a,",esquema,".itt_i_algoritmo b,",esquema,".itt_i_maestros c where b.kit=a.itm and b.itm=",    strInt(Y[2])," and a.drqj=",strInt(fecha)," and a.itm=c.itm and a.uorg>0"]),conn)
    if size(dem,1)>0
        dem=convert(Matrix{Any},dem)
        dm=convert(Matrix{Float64},dem[:,3:5])
        modele=Model(solver=CbcSolver(seconds=tiempo,ratioGap=gap,threads=tareas))
        @variable(modele,0<=x[1:size(dm,1)]<=10^5,category=:Int)
        @variable(modele,0<=e<=10^5)
        @variable(modele,0<=e1<=10^5)
        v=reshape(dm[:,3],1,size(dm,1))
        @objective(modele,Min,sum(mean(dm[:,1:2],dims=2)'.^-1*x))
        @constraint(modele,x.>=dm[:,1])
        @constraint(modele,x.-e.<=dm[:,2])
        @constraint(modele,v*x.+e1.>=Y[3])
        @constraint(modele,v*x.-e1.<=Y[3])
        @constraint(modele,sum(x)>=.9*sum(dm[:,1:2])/2)
        @constraint(modele,sum(x)-e<=1.08*sum(dm[:,1:2])/2)
        statuse=solve(modele)
        if statuse==:Optimal
            x_value=getvalue(x)
            e_value=getvalue(e)
            dem=hcat(repeat(hcat(fecha,reshape(Y,1,prod(size(Y)))),size(dem,1),1),dem,x_value)
        else
            dem=hcat(repeat(hcat(fecha,reshape(Y,1,prod(size(Y)))),size(dem,1),1),dem,repeat(["Problema Infactible de producto terminado"],size(dem,1),1))
            println("Problema Infactible de producto terminado\n")
        end
    else
        dem=hcat(fecha,reshape(Y,1,prod(size(Y))),-ones(1,5),"No hay demanda de producto terminado")
        println("No hay demanda de producto terminado\n")
    end
    return dem
end
function clearTable(tabla::String,conn::Array{String})
    dsn,username,password=conn
    db=ODBC.DSN(dsn,username,password)
    stmt1=join(["delete from ",tabla])
    qry=ODBC.query(db,join(["SELECT * FROM ",tabla]))
    if size(qry,1)>0
        ODBC.execute!(db,stmt1)
    end
    ODBC.disconnect!(db)
end
function main()
    dsn="JDE_ODBC_ORACLE"
    username="CONS_CARNICOS"
    password="C4rn1c0s"
    conn=[dsn,username,password]
    esquema="INTEGRA"
    clearTable(join([esquema,".ITT_O_MPS"]),conn)
    clearTable(join([esquema,".ITT_O_WO"]),conn)
    option2=parse(Bool,EjecutarQuery(join(["SELECT INEV01 FROM ",esquema,".ITT_I_parametros WHERE INUKID=2"]),conn).INEV01[1])
    #parametros de optimizacion
    tiempo=float(EjecutarQuery(join(["SELECT INAA/100 as num FROM ",esquema,".ITT_I_parametros WHERE INUKID=1"]),conn).NUM[1]) #segundos
    gap=float(EjecutarQuery(join(["SELECT INAA/100 as num FROM ",esquema,".ITT_I_parametros WHERE INUKID=6"]),conn).NUM[1]) #taza de aproximacion al optimo
    tareas=Int(EjecutarQuery(join(["SELECT INAA/100 as num FROM ",esquema,".ITT_I_parametros WHERE INUKID=7"]),conn).NUM[1]) #trabajos totales en paralelo
    OptimoCarnes(esquema,option2,tiempo,gap,tareas,conn)
    println("proceso finalizado\n")
    db=ODBC.DSN(dsn,username,password)
    stmt1=join(["update ",esquema,".ITT_I_Control set alev01='3'"])
    ODBC.execute!(db,stmt1)
    ODBC.disconnect!(db)
end
main()
