using Distributions, Plots, StatPlots, StatsBase
gr()

# conditional posterior function we want draws from:
function condpost_theta(theta)
    ptheta = 0.8*exp(-0.5*theta^2) + 0.1*exp((-(theta - 3.0)^2)/8)
    return ptheta
end

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
theta0 = 5.0    # starting value
M = 100000  # number of draws
sigthc = 2.0  # turning parameter

#function MHsampler(f, th_start,tuning; M=10000)
thdraw = zeros(M)
thdraw[1] = theta0

for i in 2:M
    thetac = thdraw[i-1] + sigthc*randn()  # new candidate
    pthc = condpost_theta(thetac)
    ptht = condpost_theta(thdraw[i-1])
    pp = pthc/ptht
    r = rand()
    if (r < pp)
        thdraw[i] = thetac
    else
        thdraw[i] = thdraw[i-1]
    end
end

rejrate(thdraw)

# timeplot of draws
p1 = plot(thdraw,linecolor=:green,label="MH draws")

# posterior density - burnin
mh_sample_th = thdraw[1001:M] # drop burn-in
p2 = plot(mh_sample_th,st=:density,label="MH density")

# plot ACF of draws
rho = autocor(mh_sample_th)
p3 = plot(rho,st=:scatter)
plot(p1,p2,p3,layout=(1,3),legend=false)
