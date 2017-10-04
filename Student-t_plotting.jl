# PLotting Student-t - analytical and from MC draws
using Distributions, Plots
gr()
b = -5.0:0.1:5.0
# Standardized centralized distributions
pbt = pdf(TDist(98),b)
pbn = pdf(Normal(0,1),b)
plot(b,pbt,linecolor=:blue,label="Student-t")
plot!(b,pbn, linecolor=:green,label="Normal")

# MC draws
using StatPlots
bdraws = rand(Normal(0,1),100000)
plot!(bdraws,st=:density,linecolor=:red,label="Normal draws")
bdrawst = rand(TDist(100),100000)
plot!(bdraws,st=:density,linecolor=:orange,label="t draws")

# now how to plot pdf analytically for Student-t when noncentral
# Using following values for mean and se:
b1 = 2.025
seb1 = 0.676
pbn = pdf(Normal(b1,seb1),b)
plot(b,pbn,linecolor=:blue,label="Normal analytical")

bdraws = rand(Normal(b1,seb1),100000)
plot!(bdraws,st=:density,linecolor=:red,label="Normal draws")
bdrawst = rand(TDist(100),100000).*seb1 + b1
plot!(bdrawst,st=:density,linecolor=:orange,label="t draws")

## using the formula for Student-t pdf:
function student_t(x,mu,sig,v)
    lc = lgamma((v+1)/2) - 0.5*log(pi*v*sig^2) - lgamma(v/2)
    k = (1 + ((x - mu)^2)/(v*sig^2))^(-(v+1)/2)
    return exp(lc)*k
end

pstb = zeros(length(b))
for i in 1:length(b)
    pstb[i] = student_t(b[i],b1,seb1,98)
end
plot!(b,pstb,linecolor=:green,label="analytical t formula",legend=:topleft)
