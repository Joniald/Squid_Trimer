cd("C:/Users/User/Desktop/Trimer_Squid_Git")

################################################################################

using DifferentialEquations
using Plots
using DynamicalSystems
import Statistics
using MAT

################################################################################

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

################################################################################

u0 = [4.0612, 2.4443, 4.0612, 1.2859, 4.8265, 2.7072, 0]
p = [0.024, 0.1369, 0.02, 0, 0.1075, 1.1]
ds = ContinuousDynamicalSystem(squid, u0, p)
tinteg =  tangent_integrator(ds, 7)

################################################################################

om = 1.233
period = (2*pi)/om
set_parameter!(ds,6,om)
periods = 10000
t_interval = periods*period
dtt=period*0.01
dataC = trajectory(ds, t_interval, dt = dtt, u0, Ttr = 10000*period)

intPeriods = 9000
plot( ((intPeriods*100):(periods*100))./100, dataC[(intPeriods*100):end-1,1], markersize=0.3, color = 1, xlabel="t", ylabel="Φ", label = "")
p1 = plot!( ((intPeriods*100):(periods*100))./100, dataC[(intPeriods*100):end-1,3], markersize=0.3, color = 2, xlabel="t",  ylabel="Φ", label = "")
#title = "This is Plotted using Plotly")
p3 = plot( dataC[(intPeriods*100):end-1,1], dataC[(intPeriods*100):end-1,2], markersize=0.3, color = 1, xlabel="φ2",  ylabel="φ1", label = "")
p5 = plot( dataC[(intPeriods*100):end-1,1], dataC[(intPeriods*100):end-1,3], markersize=0.3, color = 1, xlabel="φ3", ylabel="φ1", label = "")

om = 1.24
period = (2*pi)/om
set_parameter!(ds,6,om)
dataCS = trajectory(ds, t_interval, dt = dtt, u0, Ttr = 10000*period)

plot( ((intPeriods*100):(periods*100))./100, dataCS[(intPeriods*100):end-1,1], markersize=0.3, color = 1, xlabel="t", label = "")
p2 = plot!( ((intPeriods*100):(periods*100))./100, dataCS[(intPeriods*100):end-1,3], markersize=0.3, color = 2, xlabel="t", label = "")
p4 = plot( dataCS[(intPeriods*100):end-1,1], dataCS[(intPeriods*100):end-1,2], markersize=0.3, color = 1, xlabel="φ2", label = "")
p6 = plot( dataCS[(intPeriods*100):end-1,1], dataCS[(intPeriods*100):end-1,3], markersize=0.3, color = 1, xlabel="φ3", label = "")

plot(p1, p2, p3, p4, p5, p6, layout = (3, 2))
savefig("Figures/figTimeSeries.png")

################################################################################

etaParC = sqrt.( (dataC[:,1]-dataC[:,3]).*(dataC[:,1]-dataC[:,3])  +  (dataC[:,4]-dataC[:,6]).*(dataC[:,4]-dataC[:,6]))
#Statistics.mean(etaParC)
etaParCS = sqrt.( (dataCS[:,1]-dataCS[:,3]).*(dataCS[:,1]-dataCS[:,3])  +  (dataCS[:,4]-dataCS[:,6]).*(dataCS[:,4]-dataCS[:,6]))

plot((1:(periods*100))./100, etaParC[1:end-1], markersize=0.3, color = 1,
 xlabel="1/T", ylabel="eta", label = "")
plot!((1:(periods*100))./100, etaParCS[1:end-1], markersize=0.3, color = 2,
 xlabel="1/T", ylabel="eta", label = "", title = "
 η for Ω = 1.24
 ")

savefig("Figures/figMeanEta.png")
################################################################################

realizations = 1:1
om_initial = 1.22
om_final = 1.25
period_final = 2*pi/om_initial
period_initial = 2*pi/om_final
period_step = 0.001 # The step of bifurcation on period.
period_range = period_initial:period_step:period_final # The range of bifurcation on period.
n_points = 10 # The nr of points that we want to keep for each value of the stroboscopic variable.
nr_period = 1000 # The nr of periods that we run the system plus nr of periods for the transient state.
per_step = 0.01 # The proportion of period that we want as an output.

pylist_lyap = zeros(length(period_range), 7)
f1_max = zeros(length(period_range), n_points)
eta_mean = zeros(length(period_range))
omega = zeros(length(period_range))

for j = 1:length(realizations)
  s = rand(6)
  u0 = [5*s[1], 5*s[2], 5*s[3], 5*s[4], 5*s[5], 5*s[6], 0]

    for (i, period) in enumerate(period_range)

      # Solving the system
      om = (2*pi)/period
      set_parameter!(ds,6,om)
      data = trajectory(ds, nr_period*period, dt=period*per_step, u0, Ttr = nr_period*period)

      # Calculating the mean value over time of eta parameter
      eta_par = sqrt.( (data[:,1]-data[:,3]).*(data[:,1]-data[:,3])  +  (data[:,4]-data[:,6]).*(data[:,4]-data[:,6]))
      eta_mean[i] = Statistics.mean(eta_par)

      # Calculating all Lyapunov exponents
      reinit!(tinteg, u0, orthonormal(7, 7))
      pylist_lyap[i,:] = lyapunovs(tinteg, nr_period*period, 0.5*period, nr_period*period)

      # Stroboscopic calculation. The maxima in each period.
      z1 = zeros(nr_period-1)
      s1 =convert(Int64,round(1.0/per_step,digits=0))
      for i=1:nr_period-1
         z1[i] = maximum( data[ i*s1+1:(i+1)*s1,1] ) # If you want stoboscopic in Φ2 variable you have to put 2 instead of 1.
      end

      # Holding only n_points for the stroboscopic map.
      f1_max[i,:]=z1[nr_period-(n_points+1):nr_period-2]

      # We have to put some noise because the system is multistable.
      omega[i]=om
      s = 0.00001*rand(4)
      u0 = [data[end,1],data[end,2]+s[1],data[end,3],data[end,4]+s[2],data[end,5]+s[3],data[end,6]+s[4],0]
      println(i,",",j,"")
  end
end

################################################################################

plot(omega, pylist_lyap[:,1], markersize=1.0, color = 1, xlabel="Ω", ylabel="Lyapunovs", label = "L1", legend=:topleft)
plot!(omega, pylist_lyap[:,2], markersize=1.0, color = 2, xlabel="Ω", ylabel="Lyapunovs", label = "L2", legend=:topleft)
plot!(omega, pylist_lyap[:,3], markersize=1.0, color = 3, xlabel="Ω", ylabel="Lyapunovs", label = "L3", legend=:topleft)
plot!(omega, pylist_lyap[:,4], markersize=1.0, color = 4, xlabel="Ω", ylabel="Lyapunovs", label = "L4", legend=:topleft)
plot!(omega, pylist_lyap[:,5], markersize=1.0, color = 5, xlabel="Ω", ylabel="Lyapunovs", label = "L5", legend=:topleft)
plot!(omega, pylist_lyap[:,6], markersize=1.0, color = 6, xlabel="Ω", ylabel="Lyapunovs", label = "L6", legend=:topleft)
p01 = plot!(omega, pylist_lyap[:,7], markersize=1.0, color = 7, xlabel="Ω", ylabel="Lyapunovs", label = "L7", legend=:topleft)

p02 = plot(omega, eta_mean[:], markersize=1.0, color = 1, xlabel="Ω", ylabel="eta_mean", label = "")

p03 = scatter(omega, f1_max, markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")

plot(p01, p02, p03, layout = (3, 1))
savefig("Figures/figStrobLyapEta.png")

################################################################################

nr_repetitions = 15 # The number of repetitions
f01_max = zeros(length(period_range), nr_repetitions, n_points)
omega01 = zeros(length(period_range))

for j=1:nr_repetitions
    s = rand(6)
    u0 = [5*s[1], 5*s[2], 5*s[3], 5*s[4], 5*s[5], 5*s[6], 0]
 for (i, period) in enumerate(period_range)
   # Parameters.
   om = (2*pi)/period
   set_parameter!(ds,6,om)
   #global nr_period, n_points, per_step

   # Runing the system.
   data = trajectory(ds, nr_period*period, dt=period*per_step, u0, Ttr = nr_period*period)

   # In order to have the maxima in each period.
   z1 = zeros(nr_period-1)
   s1 =convert(Int64,round(1.0/per_step,digits=0))
   for i=1:nr_period-1
      z1[i] = maximum( data[ i*s1+1:(i+1)*s1,1] )
   end

   # Holding only n_points.
   f01_max[i,j,:]=z1[nr_period-(n_points+1):nr_period-2]
   omega01[i]=om

   # We have to put some noise because the system is multistable.
   s = 4*rand(4)
   u0 = [data[end,1],data[end,2]+s[1],data[end,3],data[end,4]+s[2],data[end,5]+s[3],data[end,6]+s[4],0]

    println(i,",",j,"")
 end
end

scatter(omega01, f01_max[:,1,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,2,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,3,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,4,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,5,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,6,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,7,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,8,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,9,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,10,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,11,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,12,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,13,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,14,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")
scatter!(omega01, f01_max[:,15,:], markersize=1.5, legend=:none, color = "black", xlabel="Ω", ylabel="Φ1max")


savefig("Figures/figStrob.png")

################################################################################

om_range = 1:1:104
sum_p = zeros(length(om_range))
omega = zeros(length(om_range))
for i = 1:length(sum_p)
 sum_p[i] = pylist_lyap[i,1]+pylist_lyap[i,2]+pylist_lyap[i,3]+pylist_lyap[i,4]+pylist_lyap[i,5]+pylist_lyap[i,6]+pylist_lyap[i,7]+3*p[1]
 omega[i] = i
end
omega
scatter(omega[:], sum_p[:], markersize=3.0, color = 1)


################################################################################
