using Distributions, Plots, StatPlots
gr()
#### CORRECT FORM TO DRAW FROM IG Distribution
# Suppose s^2 = 0.3, and n = 100:
a = 100/2
b = 0.3*100/2
M = 100000
s2draws = 1./rand(Gamma(a,1/b),M)
plot(s2draws, st=:density,color="blue",label="Gamma")
mean(s2draws)
std(s2draws)
median(s2draws)

# exponential with parameter Θ = 1/b
Θ = 2.0
y = rand(Exponential(1/Θ),1000)
plot(y, st=:histogram,normed = true, color=:lightblue, label="Exponential")
plot!(y, st=:density, linecolor=:orange, linewidth = 3, label="density")
mle = 1/mean(y)
