using ModelingToolkit
using NonlinearSolve


 abstract type SymbolicKineticModel <: KineticModel end

 #TODO: This sucks
 function Base.show(io::IO, m::SymbolicKineticModel)
    s = show(io,m.model)
    print(io, s)
end

 struct SymbolicButlerVolmer <: SymbolicKineticModel
    model::NonlinearSystem
 end

"""

"""
 get_model(m::SymbolicKineticModel) = m.model


"""
    SymbolicButlerVolmer(;name)

Generates a symbolic (ModelingToolkit) model of the butler volmer equation. 
"""
function SymbolicButlerVolmer(;name)
    @variables V,kT,rate
    @parameters A,α
    eqs = [rate~A*(exp(α*V/kT)-exp(-(1-α)*V/kT))]
    return SymbolicButlerVolmer(NonlinearSystem(eqs,[V,kT,rate],[A,α],name=name))
end

"""
    SymbolicButlerVolmer(bv::ButlerVolmer;name)

Generates a symbolic (ModelingToolkit) model of the butler volmer equation, with default parameter values given in the standard ElectrochemicalKinetics butler-volmer model bv.
"""
function SymbolicButlerVolmer(bv::ButlerVolmer;name)
    sys = SymbolicButlerVolmer(name=name)
    sys.model.α = bv.α
    sys.model.A = bv.A
    return sys
end

function (model::SymbolicButlerVolmer)(;V_app,kT_app)
    params = parameters(model.model)
    defaults = ModelingToolkit.defaults(model.model)
    for param in params
        try
            defaults[param]
        catch
            throw(ErrorException("Parameters need to be set"))
        end
    end
    eq = [
        model.model.V~V_app
        model.model.kT~kT_app
    ]
    ns = compose(NonlinearSystem(eq,[],[],name=:ns),model.model)
    
    guess = [
        model.model.rate=>14.0
        model.model.V=>V_app
        model.model.kT=>kT_app
    ]
    prob = NonlinearProblem(ns,guess)
    sol = solve(prob,NewtonRaphson())
end

