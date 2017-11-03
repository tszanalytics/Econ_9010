using Distributions, Plots, StatPlots, StatsBase
gr()

# conditional posterior function we want draws from:
function ln_condpost(alpha,beta,sigma,y,X)
    n = length(y)
    b = [alpha; beta]
    ptheta = -(n+1)*log(sigma^2) - sigma^(-2)*(y - X*b)'*(y - X*b)
    return ptheta
end

function RWMHsampler1(f, theta,tuning)
  thdraw = zeros(M)
  thdraw[1] = th_start
  for i in 2:M
      thetac = thdraw[i-1] + tuning*randn()  # new candidate
      pthc = f(thetac)
      ptht = f(thdraw[i-1])
      pp = pthc/ptht
      r = rand()
      if (r < pp)
        thdraw[i] = thetac[1]
      else
        thdraw[i] = thdraw[i-1]
      end
  end

M = 100000
tuninga = 2.0

adraw = zeros(M)

beta0 = 2.0
sigma = 1.0
adraw[1] = 1.0

for i in 2:M
# MH for alpha
  adraw[i] =  adraw[i-1] + tuninga*randn()
  pthc = ln_condpost(adraw[i],beta0,sigma,y,X)
  ptht = ln_condpost(adraw[i-1],beta0,sigma,y,X)
  pp = exp(pthc - ptht)
  r = rand()
  if (r < pp)
    adraw[i] = adraw[i]
  else
    adraw[i] = adraw[i-1]
  end
end

plot(adraw)
plot(adraw,st=:density)



x = randn(10)
X = [ones(10) x]
y = X*[1.0; 2.0] + randn(10)
plot(x,y, st=:scatter)

a = -1.0:0.1:3.0
m = length(a)
lpa = zeros(m)
for i in 1:m
  lpa[i] = ln_condpost_alpha(a[i],2.0,1.0,y,X)
end

pa = exp.(lpa)
plot(a,pa)


# compute rejection rate
function rejrate(mh_sample)
  m = length(mh_sample)
  rej = 0
  for i in 2:m
      if(mh_sample[i] == mh_sample[i-1])
          rej += 1
      end
  end
  return rrate = rej/m
end

# random walk (RW)metropolis (MH) algorithm
# f = cond. posterior function to evaluate
# th_start = starting value
# tuning = SD of random walk to tune step size
function RWMHsampler(f, th_start,tuning; M=10000)
  thdraw = zeros(M)
  thdraw[1] = th_start
  for i in 2:M
      thetac = thdraw[i-1] + tuning*randn()  # new candidate
      pthc = f(thetac)
      ptht = f(thdraw[i-1])
      pp = pthc/ptht
      r = rand()
      if (r < pp)
        thdraw[i] = thetac[1]
      else
        thdraw[i] = thdraw[i-1]
      end
  end
  return thdraw
end

# cmean = candidate Normal distribution mean
# cstd = candidate Normal distribution SD
# theta0 = starting value
function IndNormalMHsampler(f,cmean,cstd, theta0; M=10000)
  thdraw = zeros(M)
  thdraw[1] = theta0
  for i in 2:M
      thetac = rand(Normal(cmean,cstd),1)  # new candidate
      pthc = f(thetac[1])
      ptht = f(thdraw[i-1])
      pp = pthc/ptht
      r = rand()
      if (r < pp)
        thdraw[i] = thetac[1]
      else
        thdraw[i] = thdraw[i-1]
      end
  end
  return thdraw
end


# obtain MH MCMC sample
theta0 = 5.0    # starting value
M = 10000  # number of draws
sigthc = 1.5  # tuning parameter
thdraw = RWMHsampler(condpost_theta,theta0,sigthc,M=M)

mh_sample_th = thdraw[1001:M] # drop burn-in
[mean(mh_sample_th), std(mh_sample_th)]

quantile(mh_sample_th,[0.005,0.025,0.5,0.0975,0.995])
rejrate(thdraw)
rejrate(mh_sample_th)

# timeplot of draws
p1 = plot(thdraw,linecolor=:green,label="MH draws")
# posterior density - burnin
p2 = plot(mh_sample_th,st=:density,label="MH density")
# plot ACF of draws
rho = autocor(mh_sample_th)
#p3 = plot(rho,st=:scatter)
p3 = bar(rho,bar_width=0.1)
plot(p1,p2,p3,layout=(1,3),legend=false)


# independent normal MH
theta0 = 0.0
thdrawi = IndNormalMHsampler(condpost_theta, 0.5, 3.0, theta0,M=M)


mh_sample_thi = thdrawi[1001:M] # drop burn-in
[mean(mh_sample_thi), std(mh_sample_thi)]

quantile(mh_sample_thi,[0.005,0.025,0.5,0.0975,0.995])
rejrate(thdrawi)
rejrate(mh_sample_thi)

# timeplot of draws
p1 = plot(thdrawi,linecolor=:green,label="MH draws")
# posterior density - burnin
p2 = plot(mh_sample_thi,st=:density,label="MH density")
# plot ACF of draws
rho = autocor(mh_sample_thi)
#p3 = plot(rho,st=:scatter)
p3 = bar(rho,bar_width=0.1)
plot(p1,p2,p3,layout=(1,3),legend=false)

p4 = plot(mh_sample_thi,st=:density,label=" Ind MH density")
plot(p2,p4,layout=(1,2),legend=false)


plot!(mh_sample_th,st=:density,label=" RW MH density")
