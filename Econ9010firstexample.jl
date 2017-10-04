# Generate n draws from a Normal:
using Distributions, Plots, StatPlots
gr()
n = 100000
x = rand(Normal(-1,2),n)
mean(x)

# Histogram and density plot of draws:
plot(x, st=:histogram,normed=true,color=:purple, bins=100, label = "histogram")
plot!(x, st=:density,color=:orange, label = "density")

## cdf, pdf, p-values:
cdf(Normal(0,1),1.86)
1 - cdf(Normal(0,1),1.96)
pdf(Normal(0,1),1.96)
pdf(TDist(5),1.96)
pdf(Gamma(1,1),1.96)

## OVERLOADING FUNCTIONS
# load marg_post_mu function, two versions
# using Distributions
""""
Draws from Normal with mean =  μ, SD = σ
"""
function marg_post_mu(μ,σ; M = 10000)
  # need to either evaluate pdf of t, or
  # draw from t and plot pseudo-sample
  z = rand(Normal(μ,σ),M)
  return z
end

"""
  Compute marginal posterior draws give posterior summary statistics
  input posterior mean = μ, SD = σ, and number of obs. = n
  M = number of draws
  Default M = 10000
"""
function marg_post_mu(μ,σ,n; M = 10000)
  v = n - 1
  vs2 = v*σ^2
  # need to either evaluate pdf of t, or
  # draw from t and plot pseudo-sample
  ts = μ + σ*rand(TDist(v),M)
  return ts
end

# Draw from t and plot density
zt = marg_post_mu(2.0,1.0,5,M=20000)
plot(z, st=:density,linecolor=:orange,linewidth=3,label="t density")

# Draw from Gaussian and plot density on same figure
zn = marg_post_mu(2.0,1.0,M=20000)
plot!(zn, st=:density,linecolor=:green,label="normal density",linewidth=3)



using Plots, StatPlots, Distributions
gr()
d1 = sample(1:6,1000)
d2 = sample(1:6,1000)
plot(d1, st=:histogram, normed=true, color=:pink,alpha=0.3)
plot!(d2, st=:histogram, normed=true, color=:blue,alpha=0.3)
plot!(d2, st=:density, normed=true, color=:blue,alpha=0.3)
plot!(d1, st=:density, normed=true, color=:blue,alpha=0.3)

# CLT in action - at its best!
M =10000
d1 = sample(1:6,M)
d2 = sample(1:6,M)
d3 = sample(1:6,M)
d4 = sample(1:6,M)
d5 = sample(1:6,M)
d6 = sample(1:6,M)
plot(d1, st=:histogram, normed=true, linecolor=:blue,linewidth=2)
plot(d1+d2, st=:histogram, normed=true, linecolor=:blue,linewidth=2)
plot(d1+d2+d3, st=:histogram, normed=true, linecolor=:blue,linewidth=2)
plot(d1+d2+d3+d4, st=:histogram, normed=true, linecolor=:blue,linewidth=2)
plot(d1+d2+d3+d4+d5, st=:histogram, normed=true, linecolor=:blue,linewidth=2)
