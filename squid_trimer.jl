cd("C:/Users/User/Desktop/Trimer_Squid_Git")


using DifferentialEquations
using Plots
using DynamicalSystems
import Statistics
using MAT


@inline @inbounds function squid(u, p, t)
   g, v, fac, fdc, e, om = p
   x1, x2, x3, y1, y2, y3, z = u
   du1 = y1
   du2 = y2
   du3 = y3
   du4 = -g*y1 -v*sin(2*pi*x1) + (64/(64-129*e*e+16*e*e*e))*( (1-e*e-e+(e*e/8)+e*e-(e/8))*(fac*cos(z) + fdc)      -(1-e*e)*x1 -(-e+(e*e/8))*x2 -(e*e-(e/8))*x3 )
   du5 = -g*y2 -v*sin(2*pi*x2) + (64/(64-129*e*e+16*e*e*e))*( (-e+(e*e/8)+1-(e*e/64)-e+(e*e/8))*(fac*cos(z) + fdc) -(-e+(e*e/8))*x1 -(1-(e*e/64))*x2 -(-e+(e*e/8))*x3 )
   du6 = -g*y3 -v*sin(2*pi*x3) + (64/(64-129*e*e+16*e*e*e))*( (e*e-(e/8)-e+(e*e/8)+1-e*e)*(fac*cos(z) + fdc)      -(e*e-(e/8))*x1 -(-e+(e*e/8))*x2 -(1-e*e)*x3 )
   du7 = om
   return SVector{7}(du1, du2, du3, du4, du5, du6, du7)
end

u0 = [0.0612, 2.4443, 0.0612, 1.2859, 3.8265, 2.7072, 0]
p = [0.024, 0.1369, 0.02, 0, 0.1075, 1.1]
ds = ContinuousDynamicalSystem(squid, u0, p)
tinteg =  tangent_integrator(ds, 7)



om = 1.24
period = (2*pi)/om
set_parameter!(ds,6,om)
t_interval = 5000*period
dtt=period*0.01
data = trajectory(ds, t_interval, dt = dtt, u0, Ttr = 10000*period)
eta_par = sqrt.( (data[:,1]-data[:,3]).*(data[:,1]-data[:,3])  +  (data[:,4]-data[:,6]).*(data[:,4]-data[:,6]))
Statistics.mean(eta_par)
(1:(100*100))./100
plot!((1:(5000*100))./100, eta_par[1:end-1], markersize=0.3, color = 2, xlabel="1/T", ylabel="eta", label = "Ω = 1.24")

plot((9001:(100*100))./100, data[9001:end-1,1], markersize=0.3, color = 1, xlabel="t", ylabel="Φ")
plot!((9001:(100*100))./100, data[9001:end-1,2], markersize=0.3, color = 2, xlabel="t", ylabel="Φ")
plot!((9001:(100*100))./100, data[9001:end-1,3], markersize=0.3, color = 3, xlabel="t", ylabel="Φ")

savefig("Figures/Fig_02.png")
#corr1 = Statistics.cor(data[2500000:end,1],data[2500000:end,3])
#corr2 = Statistics.cor(data[2500000:end,1],data[2500000:end,2])
#corr3 = Statistics.cor(data[2500000:end,2],data[2500000:end,3])
#plot(data[2500000:3500000,2],data[2500000:3500000,5])
#plot!(data[2500000:3500000,3],data[2500000:3500000,6])
#plot!(data[2500000:3500000,1],data[2500000:3500000,4])
#plot!(data[3000000:end,1])
#plot(data[3000000:end,3])
#plot!(data[2500000:end,2])

#z1=data[1:end,1]
#z3=data[1:end,3]
#z4=data[1:end,4]
#z6=data[1:end,6]
#z2=data[1:end,2]
#file = matopen("datat1.mat", "w")
#write(file, "datat1", z1)
#close(file)
#file = matopen("datat3.mat", "w")
#write(file, "datat3", z3)
#close(file)
#file = matopen("datat4.mat", "w")
#write(file, "datat4", z4)
#close(file)
#file = matopen("datat6.mat", "w")
#write(file, "datat6", z6)
#close(file)
#file = matopen("datat2.mat", "w")
#write(file, "datat2", z2)
#close(file)


#reinit!(tinteg, u0, orthonormal(7, 7))
#tinteg =  tangent_integrator(ds, 7)
#lyap=lyapunovs(tinteg, 3500*period, 1*period, 1500*period)
#reinit!(tinteg, u0, orthonormal(7, 7))
#gal = gali(tinteg, 5000, 5, 1e-15)[2][end]


realizations = 1:1
om = 1.245
period = (2*pi)/om
#per=5.15015:-0.0001:5.04673;
per=5.04673:0.001:5.15015;
pylist1 = zeros(length(per), 100)
omega = zeros(length(per))
pylist_lyap = zeros(length(per), 7)
#pylist2 = zeros(length(per), 100)
#gal = zeros(length(per));
eta_mean = zeros(length(per))
#initial_cond = zeros(length(per), 7)


for j = 1:length(realizations)
  s = rand(6)
  u0 = [5*s[1], 5*s[2], 5*s[3], 5*s[4], 5*s[5], 5*s[6], 0]

    for (i, period) in enumerate(per)

      om = (2*pi)/period
      set_parameter!(ds,6,om)
      #initial_cond[i,:] = u0
      data = trajectory(ds, 1000*period, dt=period*0.01, u0, Ttr = 1000*period)

      eta_par = sqrt.( (data[:,1]-data[:,3]).*(data[:,1]-data[:,3])  +  (data[:,4]-data[:,6]).*(data[:,4]-data[:,6]))
      eta_mean[i] = Statistics.mean(eta_par)

      #corr[i] = Statistics.cor(data[:,1],data[:,3])

      reinit!(tinteg, u0, orthonormal(7, 7))
      pylist_lyap[i,:] = lyapunovs(tinteg, 1000*period, 0.5*period, 1000*period)


      #reinit!(tinteg, u0, orthonormal(7, 7))
      #gal[i] = gali(tinteg, 5000, 5, 1e-15)[2][end]

      #z1=data[1:1000:end,1]
      #pylist1[i,:]=z1[3001:3100]
      omega[i]=om
      #u0 = [data[end,1],data[end,2],data[end,3],data[end,4],0]
      s = 0.00001*rand(4)
      u0 = [data[end,1],data[end,2]+s[1],data[end,3],data[end,4]+s[2],data[end,5]+s[3],data[end,6]+s[4],0]

      println(i,",",j,"")
  end
end

plot(omega, pylist_lyap[:,1], markersize=1.0, color = 1, xlabel="Ω", ylabel="Lyapunovs", label = "L1", legend=:topleft)
plot!(omega, pylist_lyap[:,2], markersize=1.0, color = 2, xlabel="Ω", ylabel="Lyapunovs", label = "L2", legend=:topleft)
plot!(omega, pylist_lyap[:,3], markersize=1.0, color = 3, xlabel="Ω", ylabel="Lyapunovs", label = "L3", legend=:topleft)
plot!(omega, pylist_lyap[:,4], markersize=1.0, color = 4, xlabel="Ω", ylabel="Lyapunovs", label = "L4", legend=:topleft)
plot!(omega, pylist_lyap[:,5], markersize=1.0, color = 5, xlabel="Ω", ylabel="Lyapunovs", label = "L5", legend=:topleft)
plot!(omega, pylist_lyap[:,6], markersize=1.0, color = 6, xlabel="Ω", ylabel="Lyapunovs", label = "L6", legend=:topleft)
plot!(omega, pylist_lyap[:,7], markersize=1.0, color = 7, xlabel="Ω", ylabel="Lyapunovs", label = "L7", legend=:topleft)

pylist_lyap[:,7]
ss[1]
savefig("Fig_01.png") #if you want to save the fig.

plot!(omega, eta_mean[:], markersize=1.0, color = 1, xlabel="Ω", ylabel="eta_mean", label = "")
savefig("Fig_03.png") #if you want to save the fig.

om_range = 1:1:104
sum_p = zeros(length(om_range))
omega = zeros(length(om_range))
for i = 1:length(sum_p)
 sum_p[i] = pylist_lyap[i,1]+pylist_lyap[i,2]+pylist_lyap[i,3]+pylist_lyap[i,4]+pylist_lyap[i,5]+pylist_lyap[i,6]+pylist_lyap[i,7]+3*p[1]
 omega[i] = i
end
omega
scatter(omega[:], sum_p[:], markersize=3.0, color = 1)


palette(:tab10)
file = matopen("omega.mat", "w")
write(file, "omega", omega)
close(file)

file = matopen("yt.mat", "w")
write(file, "yt", pylist1)
close(file)

file = matopen("galt.mat", "w")
write(file, "galt", gal)
close(file)

file = matopen("lyapt.mat", "w")
write(file, "lyapt", pylist_lyap)
close(file)

file = matopen("corrt.mat", "w")
write(file, "corrt", corr)
close(file)

file = matopen("ini_cont.mat", "w")
write(file, "ini_cont", initial_cond)
close(file)
