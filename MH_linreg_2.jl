using Distributions, Plots, StatPlots, StatsBase
gr()


srand(1729)
x = rand(Normal(2,2),100)
e = randn(100)
y = 2 + 3*x + e

plot(x,y,st=:scatter)

X = [ones(length(y)) x]

function ln_condpost(alp,bet,sig,y,X)
  n = length(y)
  b = [alp; bet]
  palpha = -(n+1)*log(sig^2) - (sig^-2)*(y-X*b)'*(y-X*b)
  return palpha
end

M = 100000
alphastart = 10.0
betastart = 5.
sigstart = 3.
tune = 0.001:0.01:0.51
m = length(tune)
alphasample = zeros(M,m)
betasample = zeros(M,m)
sigmasample = zeros(M,m)
alphadraw = zeros(M)
betadraw = zeros(M)
sigdraw = zeros(M)
alphadraw[1] = alphastart
betadraw[1]= betastart
sigdraw[1] = sigstart
for j in 1:m
for i in 2:M
    alphac = alphadraw[i-1] + tune[j]*randn()  # new candidate
    pthc = ln_condpost(alphac,betadraw[i-1],sigdraw[i-1],y,X)
    ptht = ln_condpost(alphadraw[i-1],betadraw[i-1],sigdraw[i-1],y,X)
    pp = exp(pthc - ptht)
    r = rand()
    if (r < pp)
      alphadraw[i] = alphac
    else
      alphadraw[i] = alphadraw[i-1]
    end
    betac = betadraw[i-1] + tune[j]*randn() # new candidate for beta
    bthc = ln_condpost(alphadraw[i],betac,sigdraw[i-1],y,X)
    btht = ln_condpost(alphadraw[i],betadraw[i-1],sigdraw[i-1],y,X)
    bb = exp(bthc - btht)
    rbeta = rand()
    if (rbeta < bb)
        betadraw[i] = betac
    else
        betadraw[i] = betadraw[i-1]
    end
    sigc = sigdraw[i-1] + tune[j]*randn() # new candidate for sigma
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
alphasample[:,j] = alphadraw
betasample[:,j] = betadraw
sigmasample[:,j] = sigdraw
end


q = length(alphasample[1001:M,1])
rejectratealpha= zeros(m)
for p in 1:m
reject=0
for i in 1001:M
  if alphasample[i,p] == alphasample[i-1,p]
    reject += 1
  end
end
rejectratealpha[p]=reject
end
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
