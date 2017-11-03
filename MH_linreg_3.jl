using Distributions, Plots, StatPlots, StatsBase
gr()


srand(1729)
N = 1000
x = rand(Normal(2,2),N)
u = randn(N)
y = 2 + 3*x + u

truyhat = 2 + 3*x
trusighat = sum((y - truyhat).^2)/N
trutauhat = 1/trusighat

plot(x,y,st=:scatter)

X = [ones(length(y)) x]

function ln_condpost(alp,bet,tau,y,X)
  n = length(y)
  b = [alp; bet]
  palpha = (n-2)*log(tau) - tau*(y-X*b)'*(y-X*b)
  return palpha
end

M = 300000
alphastart = 1.0
betastart = 1.
sigstart = 1.
tunea = 0.05
tuneb = 0.05
tunes = 0.1

alphadraw = zeros(M)
betadraw = zeros(M)
sigdraw = zeros(M)
alphadraw[1] = alphastart
betadraw[1]= betastart
sigdraw[1] = sigstart

# MCMC loop
for i in 2:M

# MH step for alpha
    alphac = alphadraw[i-1] + tunea*randn()  # new candidate
    pthc = ln_condpost(alphac,betadraw[i-1],sigdraw[i-1],y,X)
    ptht = ln_condpost(alphadraw[i-1],betadraw[i-1],sigdraw[i-1],y,X)
    pp = exp(pthc - ptht)
    r = rand()
    if (r < pp)
      alphadraw[i] = alphac
    else
      alphadraw[i] = alphadraw[i-1]
    end

# MH step for beta
    betac = betadraw[i-1] + tuneb*randn() # new candidate for beta
    bthc = ln_condpost(alphadraw[i],betac,sigdraw[i-1],y,X)
    btht = ln_condpost(alphadraw[i],betadraw[i-1],sigdraw[i-1],y,X)
    bb = exp(bthc - btht)
    rbeta = rand()
    if (rbeta < bb)
        betadraw[i] = betac
    else
        betadraw[i] = betadraw[i-1]
    end

# MH step for sigma
    sigc = sigdraw[i-1] + tunes*randn() # new candidate for sigma
    if sigc <= 0.0
      sigc = abs(sigc)
    end
    sthc = ln_condpost(alphadraw[i],betadraw[i],sigc,y,X)
    stht = ln_condpost(alphadraw[i],betadraw[i],sigdraw[i-1],y,X)
    ss = exp(sthc - stht)
    rsig = rand()
    if (rsig < ss)
        sigdraw[i] = sigc
    else
        sigdraw[i] = sigdraw[i-1]
    end
end





function rejrate(mcdraw;burnin=1000)
  M = length(mcdraw)
  reject=0.
  for i in (burnin+1):M
    if mcdraw[i] == mcdraw[i-1]
      reject += 1.
    end
  end
  return reject/(M-1000)
end


b = 1001 # burnin
p1 = plot(alphadraw[b:end],st=:density,label="alpha")
p2 = plot(betadraw[b:end],st=:density,label="beta")
p3 = plot(sigdraw[b:end],st=:density,label="tau")
p4 = plot(1./sigdraw[b:end],st=:density,label="sigma")
plot(p1,p2,p3,p4,overlay=(1,3))
vline!([2.0 3.0 trutauhat trusighat],label="")
mean(alphadraw[b:end])
mean(betadraw[b:end])
mean(sigdraw[b:end])

rejrate(alphadraw)
rejrate(betadraw)
rejrate(sigdraw)

function acf_plot(mcdraw)
  rho = autocor(mcdraw)
  bar(rho,bar_width=0.1)
end

p1 = acf_plot(alphadraw[b:end])
p2 = acf_plot(betadraw[b:end])
p3 = acf_plot(sigdraw[b:end])
plot(p1,p2,p3,overlay=(1,3),label="")


p1 = plot(alphadraw[b:end],label="alpha")
p2 = plot(betadraw[b:end],label="beta")
p3 = plot(sigdraw[b:end],label="tau")
plot(p1,p2,p3,overlay=(1,3),label="")
hline!([2.0 3.0 trutauhat],label="")

rratealpha =rejectratealpha./q

q = length(alphasample[1001:M,1])
rejectratebeta= zeros(m)
for p in 1:m
reject=0
for i in 1001:M
  if betasample[i,p] == betasample[i-1,p]
    reject += 1
  end
end
rejectratebeta[p]=reject
end
rratebeta =rejectratebeta./q

q = length(sigmasample[1001:M,1])
rejectratesig = zeros(m)
for p in 1:m
reject=0
for i in 1001:M
  if sigmasample[i,p] == sigmasample[i-1,p]
    reject += 1
  end

rejectratesig[p]=reject
end
rratesig =rejectratesig./q

plot(tune,rratesig,st=:line,xlims=(0,0.55),xticks=0:0.05:0.55,ylims=(0,1),yticks=0:0.1:1,label="Rejection Rate for Sigma",ylabel="Rejection Rate",xlabel="Tuning Parameter")
plot!(tune,rratealpha,color="red",st=:line,xlims=(0,0.55),xticks=0:0.05:0.55,ylims=(0,1),label="Rejection Rate for Alpha",yticks=0:0.1:1,ylabel="Rejection Rate",xlabel="Tuning Parameter")
plot!(tune,rratebeta,color="green",st=:line,xlims=(0,0.55),xticks=0:0.05:0.55,ylims=(0,1),label="Rejection Rate for Beta",yticks=0:0.1:1,ylabel="Rejection Rate",xlabel="Tuning Parameter")
hline!([0.3],color="orange")


mean(alphadraw)
plot(alphadraw[1000:end],st=:density)

mean(betadraw)
mean(sigdraw)

rho = autocor(betadraw[1001:end])
bar(rho,bar_width=0.3)
