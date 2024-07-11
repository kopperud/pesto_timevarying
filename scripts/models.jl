using Revise
using Pesto
using CairoMakie


## read the primates tree
phy = readtree("data/primates.tre")
ρ = 0.635
primates = SSEdata(phy, ρ)

## create a silly model
λ1 = [0.3, 0.1]
μ1 = [0.15, 0.035]
η1 = 0.2 / 50

model1 = SSEconstant(λ1, μ1, η1)
rates = birth_death_shift(model1, primates)

logl1 = logL_root(model1, primates)

## plot the netdiv on the tree
treeplot(primates, rates)

## create another model
λ2 = t -> [
    0.15,
    0.2 - (t*0.002),
]
μ2 = t -> [
    0.1,
    0.2 - (t*0.002),
]
η2 = t -> 0.2 / 50

x = collect(range(0.0, maximum(phy.branching_times), length = 100))

μ2.(x)
λ2.(x)

model2 = Pesto.SSEtimevarying(λ2, μ2, η2)
rates = birth_death_shift(model2, primates)
logl2 = logL_root(model2, primates)
treeplot(primates, rates)

## non-identifiability
##
## two models (m1 and m2) are not idenfitiable from each other, **if**
## their pulled speciation rates are equal (λp1 = λp2)
##
## λp(t) = λ(t) ⊙ (1 - E(t)),
##
## where
##
## dE/dt = differential equation from BiSSE etc.


## let's calculate λp(t)

function pulled_speciation(model::SSEconstant, data)
    λ = model.λ
    ## first E(t)
    E = Pesto.extinction_probability(model, data);

    λp = t -> λ .* ( 1 .- E(t) )
    return (λp)
end

function pulled_speciation(model::Pesto.SSEtimevarying, data)
    λ = model.λ
    ## first E(t)
    E = Pesto.extinction_probability(model, data);

    λp = t -> λ(t) .* ( 1 .- E(t) )
    return (λp)
end



λ1p = pulled_speciation(model1, primates)
λ2p = pulled_speciation(model2, primates)


using StatsPlots

p = StatsPlots.plot(xflip = true, xlabel = "time before present", ylabel = "pulled sp")
y = transpose(hcat(λ1p.(x)...))
StatsPlots.plot!(p, x, y, label = "model 1")
y = transpose(hcat(λ2p.(x)...))
StatsPlots.plot!(p, x, y, label = "model 2")



## λ1p and λ2p are not equal. How can we set up a model `model3` 
## such that λ1p and λp3 are equal?
## suppose we set up λ3(t) and μ3(t) are piecewise linear functions,
## with say 50 knots over the time span of the phylogeny

using DataInterpolations
n_knots = 50
knots = collect(range(0.0, maximum(phy.branching_times), length = n_knots))
ys = rand(50)

f = DataInterpolations.LinearInterpolation(ys, knots; extrapolate = true)

f(65)
StatsPlots.plot(knots, f.(knots))

## the idea is to use an optimization algorithm to find the `ys`
## that minimizes the error in λp between the two models

function mean_squared_error(
    model1::SSEconstant,
    model2::Pesto.SSEtimevarying,
    data::SSEdata,
    knots::Vector{Float64},
    )
    n_knots = length(knots)

    λ1p = pulled_speciation(model1, data)
    λ2p = pulled_speciation(model2, data)

    deviation = vcat(λ1p.(knots)...) .- vcat(λ2p.(knots)...)
    mse = sum(deviation .^2) / n_knots
    return(mse)
end


mean_squared_error(model1, model2, primates, knots)


vcat(λ1p.(knots)...)















