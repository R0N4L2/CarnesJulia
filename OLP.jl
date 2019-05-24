using Cbc
#using GLPKMathProgInterface
#using Clp
#using AmplNLWriter
using Ipopt
using JuMP
using Clp
using Juniper
#funcion de optimizacion lineal
function LP(F::Array{<:Float64,1},Ad::Array{<:Float64,2},g::Array{<:Float64,1},c::Array{<:Float64,1},d::Array{<:Float64,1},R0::Array{<:Float64,2},R::Array{<:Float64,2},Q::Array{<:Float64,2},lam::Float64,tiprest::String)
    #F: vector de funcion objetivo
    #Ad: Matriz de restriccion de demanda
    #g: vector de demanda a cumplir
    #c: vector de cota inferior de nivel de servicio
    #b: vector de cota superoro de sobreoferta por producto
    #R: Matriz de equilibrio de peso
    #R0: Matriz de restricciones de demanda 0
    #Q: Matriz de diferencias entre priernas y brazos, para compra por cadenas
    #lam: peso que se le da a la restirccion de compra por cadenas
    #tiprest: ve si la matriz de demanda 0 es desigualdad o igualdad
    model = Model(solver =ClpSolver()) #GLPKSolverLP()
    @variable(model,0<=x[1:size(R0,2)]<=10000)
    @variable(model,0<=e[1:size(Q,1)]<=10000)
    @objective(model, Max,F'*x-lam*sum(e))
    MI=vcat(Ad,(sum(Ad,dims=1)+sum(R0,dims=1)))
    BI=vcat(g.*c,.9 .*sum(g))
    BS=vcat(g.*d,1.08.*sum(g))
    #BS=vcat(g.*b,1.08 .*sum(g))
    @constraint(model,MI*x.>=BI)
    @constraint(model,MI*x.<=BS)
    if tiprest=="igualdad"
        @constraint(model,R0*x.==0)
    else
        @constraint(model,R0*x.>=0)
    end
    @constraint(model,R*x.==0)
    @constraint(model,Q*x-e.<=0)
    @constraint(model,Q*x+e.>=0)
    status=solve(model)
    x_value=getvalue(x)
    e_value=getvalue(e)
    return status,x_value,e_value
end
#funcion de optimizacion entera
function IP(Rel::BitArray{1},F::Array{<:Float64,1},Ad::Array{<:Float64,2},g::Array{<:Float64,1},c::Array{<:Float64,1},R0::Array{<:Float64,2},R::Array{<:Float64,2},Q::Array{<:Float64,2},lam::Float64,tiprest::String)
    #Rel: son las variables relajadas osea continuas, las otras seran enteras.
    #F: vector de funcion objetivo
    #Ad: Matriz de restriccion de demanda
    #g: vector de demanda a cumplir
    #c: vector de cota inferior de nivel de servicio
    #b: vector de cota superoro de sobreoferta por producto
    #R: Matriz de equilibrio de peso
    #R0: Matriz de restricciones de demanda 0
    #Q: Matriz de diferencias entre priernas y brazos, para compra por cadenas
    #lam: peso que se le da a la restirccion de compra por cadenas
    #tiprest: ve si la matriz de demanda 0 es desigualdad o igualdad
    modele = Model(solver=CbcSolver())
    @variable(modele,0<=ve[1:sum(.!Rel)]<=10000,category=:Int)
    @variable(modele,0<=vr[1:sum(Rel)]<=10000)
    @variable(modele,0<=e[1:size(Q,1)]<=10000)
    MEe=R[:,.!Rel]
    MEr=R[:,Rel]
    Qe=Q[:,.!Rel]
    MEQ0r=R0[:,Rel]
    MEQ0e=R0[:,.!Rel]
    Wr=W[:,Rel]
    We=W[:,.!Rel]
    fe=F[.!Rel]
    fr=F[Rel]
    @objective(modele, Max,fr'*vr+fe'*ve-lam*sum(e))
    MI=vcat(Ad,(sum(Ad,dims=1)+sum(R0,dims=1)))
    MIe=MI[:,.!Rel]
    MIr=MI[:,Rel]
    BI=vcat(g.*c,.9 .*sum(g))
    BS=1.08.*vcat(g,sum(g))
    #BS=vcat(g.*b,1.08*sum(g))
    @constraint(modele,MIr*vr+MIe*ve.>=BI)
    @constraint(modele,MIr*vr+MIe*ve.<=BS)
    if tiprest=="igualdad"
        @constraint(modele,MEQ0r*vr+MEQ0e*ve.==0)
    else
        @constraint(modele,MEQ0r*vr+MEQ0e*ve.>=0)
    end
    @constraint(modele,MEr*vr+MEe*ve.==0)
    @constraint(modele,Qe*ve-e.<=0)
    @constraint(modele,Qe*ve+e.>=0)
    statuse=solve(modele)
    vr_value=getvalue(vr)
    ve_value=getvalue(ve)
    e_value=getvalue(e)
    X=ones(length(Rel))
    if statuse==:Optimal
        X[Rel]=vr_value
        X[.!Rel]=ve_value
    else
        X=NaN*X
    end
    return statuse,X,e_value
end
