#!/usr/bin/env python3
import time, math
import matplotlib.pyplot as plt
import numpy
from matplotlib.patches import Circle

def getDensity(alt):
	# Returns International Standard Atmosphere 1976 density as a function
	# of height in meters
	# https://en.wikipedia.org/wiki/Barometric_formula#Density_equations
	R_star = 8.3144598
	g0     = 9.80665
	M      = 0.0289644
	if alt < 11000:
		hb   = 0
		rhob = 1.225
		Tb   = 288.15
		Lb   = -0.0065
		rho = rhob * ( Tb / ( Tb + Lb*(alt - hb) ) )**( 1 + ( g0 * M / (R_star * Lb)))
	elif alt < 20000:
		hb   = 11000
		rhob = 0.36391
		Tb   = 216.65
		rho = rhob * math.exp( -g0 * M * (alt - hb)/(R_star * Tb) )
	elif alt < 32000:
		hb   = 20000
		rhob = 0.08803
		Tb   = 216.65
		Lb   = 0.001
		rho = rhob * ( Tb / ( Tb + Lb*(alt-hb) ) )**( 1 + ( g0 * M / (R_star * Lb)))
	elif alt < 47000:
		hb   = 32000
		rhob = 0.01322
		Tb   = 228.65
		Lb   = 0.0028
		rho = rhob * ( Tb / ( Tb + Lb*(alt-hb) ) )**( 1 + ( g0 * M / (R_star * Lb)))
	elif alt < 51000:
		hb   = 47000
		rhob = 0.00143
		Tb   = 270.65
		rho = rhob * math.exp( -g0 * M * (alt - hb)/(R_star * Tb) )
	elif alt < 71000:
		hb   = 51000
		rhob = 0.00086
		Tb   = 270.65
		Lb   = -0.0028
		rho = rhob * ( Tb / ( Tb + Lb*(alt-hb) ) )**( 1 + ( g0 * M / (R_star * Lb)))
	elif alt < 84852:
		hb   = 71000
		rhob = 0.000064
		Tb   = 214.65
		Lb   = -0.002
		rho = rhob * ( Tb / ( Tb + Lb*(alt-hb) ) )**( 1 + ( g0 * M / (R_star * Lb)))
	else:
		rho = 0
	return rho

def getMass(t):
	if t < t_meco:
		# 1st stage burn
		m = m1+m2+m3 + mp1*(t_meco-t)/tb1 + mp2
	elif t >= t_meco and t < t_sep1:
		# After MECO before SEP 1
		m = m1+m2+m3+mp2
	elif t >= t_sep1 and t < t_ssi:
		# After SEP 1 before SSI
		m = m2+m3+mp2
	elif t >= t_ssi and t < t_seco:
		# 2nd stage burn
		m = m2+m3 + mp2*(t_seco-t)/tb2
	elif t >= t_seco and t < t_sep2:
		# After SECO before SEP 2
		m = m2+m3
	else:
		m = m3
	return m

F=0

def getThrust(F_old,qmax,q,Kp,Ki,Kd,integral):
	if t<t_meco: # First stage
		# 1. Find needed change
		# 2. Set change to max if too big
		# 3. Set thrust to thrust plus change
		
		error = q-qmax
		integral = integral + error*dt
		derivative = (error - error_old)*dt
		throttle = 100-min(100,max(0,Kp*error + Ki*integral + Kd*derivative))
		
		dF = F1 * throttle/100 - F_old
		
		if dF > dFdt * dt:
			dF = dFdt * dt
		
		F = F_old + dF
		
		# list_pid_p.append(Kp*error)
		# list_pid_i.append(Ki*integral)
		# list_pid_d.append(Kd*derivative)
	elif t>t_ssi and t<t_seco: # Second stage
		F=F2
	else:
		F=0
	return F

# Environment

G  = 6.67e-11 # Gravitational constant
M  = 5.97e+24  # Earth mass in kg
mu = G*M
R0 = 6378000
g0 = mu/R0**2

# Falcon 9 Full Thrust Block 5
# https://en.wikipedia.org/wiki/Falcon_9_Full_Thrust#Block_5

m1 = 22200 # Mass of empty first stage
m2 =  4000 # Mass of empty second stage
m3 = 18000 # Mass of payload

mp1 = 433100-m1 # Mass of propellant, first stage
mp2 = 111500-m2 # Mass of propellant, second stage

tb1 = 162 # First stage burn time
tb2 = 397 # Second stage burn time

tbi1 = 5 # Time between MECO and 1st stage separation
tbi2 = 5 # Time between 1st stage separation and 2nd stage ignition
tbi3 = 60 # Time between SECO and 2nd stage separation

t_meco = tb1
t_sep1 = tb1 + tbi1
t_ssi  = tb1 + tbi1 + tbi2
t_seco = tb1 + tbi1 + tbi2 + tb2
t_sep2 = tb1 + tbi1 + tbi2 + tb2 + tbi3

# Thrust control

Kp=1
Ki=0
Kd=200

qmax = 25000 # Maximum dynamic pressure

F1 = 7607000
F2 = 934000
dFdt = 2000000 # Maximum thrust change rate [N/s]

# Pitch maneuver
t_pitch_begin    = 16 # Time for pitch maneuver begin
t_pitch_duration = 4.5 # Duration of pitch maneuver
pitch_angle      = 1 # Nozzle deflection in degrees during pitch maneuver

dia = 3.66	# Rocket diameter [m]
cd = 0.65  # Rocket drag coefficient [·]

t  = 0
dt = 0.01

m = m1 + m2 + m3 + mp1 + mp2 # Total mass

liftoff = False
maxq1   = False
maxq2   = False
meco    = False
sep1    = False
ssi     = False # Second stage ignition
seco    = False
sep2    = False
apogee  = False
perigee = False
impact  = False
orbit   = False
oncearound = False

A = math.pi * (dia/2)**2

tmax = 3*3600

# Important note about coordinate system:
# We are looking down on the earth from above the north pole.
# x is to the right, and y is up
# Rocket is launched eastward, which is counter-clockwise
# (positive direction of rotation)

#s=0 # southing
e=0 # easting NOT IN USE
z=0 # zenith  NOT IN USE

x=R0
y=0

r=R0 # distance to earth center

alt=0 # height
d=0 # downrange

vx=0
vy=0
v=0

# How to add the rotation of the earth (465 m/s at equator) but
# NOT to the velocity, so as to induce an enormous drag?

theta = 0 # position vector angle
phi   = 0 # velocity vector angle
theta_old = 0

integral=0
error=0
error_old=0

list_x = []
list_y = []
list_v = []
list_q = []
list_d = []
list_t = []
list_D = []
list_m = []
list_F = []
list_rp = []
list_alt = []
list_rho = []
list_fpa = []
list_event = []
list_event_x = []
list_event_y = []
list_event_d = []
list_event_a = []

list_pid_p = []
list_pid_i = []
list_pid_d = []

list_pid_p.append(0)
list_pid_i.append(0)
list_pid_d.append(0)

while not (impact or oncearound):

	m = getMass(t)
	
	# Drag
	rho = getDensity(alt)
	q = 0.5 * rho * v**2
	D = q * A * cd

	# Local gravitational acceleration
	g = mu/r**2

	# Angle of position vector
	theta_old = theta
	theta = math.atan2(y,x)

	# Angle of velocity vector
	phi   = math.atan2(vy,vx)

	# Flight-path-angle
	if v!=0:
		fpa = math.asin((x*vx+y*vy)/(r*v))
	else:
		fpa=0

	# Pitch maneuver
	if t>t_pitch_begin and t<t_pitch_begin+t_pitch_duration:
		phi=phi+pitch_angle/180*math.pi

	# Determine thrust
	F_old = F
	F = getThrust(F_old,qmax,q,Kp,Ki,Kd,integral)

	# Calculate the Force, Luke
	Fx = F*math.cos(phi) - m*g*math.cos(theta) - D*math.cos(phi)
	Fy = F*math.sin(phi) - m*g*math.sin(theta) - D*math.sin(phi)

	# Equations of motion
	ax = Fx/m
	ay = Fy/m

	a = math.sqrt(ax**2+ay**2)

	vx = vx+ax*dt
	vy = vy+ay*dt

	v = math.sqrt(vx**2+vy**2)

	x = x + vx*dt
	y = y + vy*dt

	r = math.sqrt(x**2+y**2)
	
	alt_old = alt
	alt = r-R0
	d = R0*theta

	sma = 1/(2/r-v**2/mu)                                  # Semimajor axis
	E = v**2/2-mu/r                                        # Specific energy
	h = r * v * math.sin( math.acos( (x*vx+y*vy)/(r*v) ) ) # Specific angular momentum
	e = math.sqrt( 1 + 2*E*h**2/mu**2 )                    # Eccentricity

	ra = sma*(1+e)
	rp = sma*(1-e)

	list_x.append(x/1000)
	list_y.append(y/1000)
	list_v.append(v)
	list_d.append(d/1000)
	list_F.append(F/1000)
	list_q.append(q)
	list_t.append(t)
	list_D.append(D)
	list_m.append(m/1000)
	list_rp.append((rp-R0)/1000)
	list_alt.append(alt/1000)
	list_rho.append(rho)
	list_fpa.append(fpa/math.pi*180)

	# Pad
	if alt<0 and not liftoff:
		x=R0
		y=0
		vx=0
		alt=0
		r=R0

	# OUTPUT
	#print("t=%.2f, alt=%.2f, d=%.2f, v=%.2f, a=%.2f, F=%.2f" % (t,alt,d,v,a,F))
	#time.sleep(0.001)
	#print("hp",round((rp-R0)/1000),"km - ha",round((ra-R0)/1000),"km")
	#print("hp %.2f - ha %.2f" % (rp-R0,ra-R0))
	#print( "%.2f, %.2f, %.2f" % (t,F,alt))
	#time.sleep(0.1)

	# Liftoff
	if alt>0 and not liftoff:
		liftoff = True
		print("Liftoff time:    ", round(t,3), "s")

	# Max Q
	if error>0 and not maxq1:
		maxq1 = True
		print("Max. Q entered:  ", round(t,3), "s")
		print("Max. Q height:   ", round(alt/1000,3), "km")
		print("Max. Q speed:    ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("Max Q enter")

	if error==0 and maxq1 and not maxq2:
		maxq2 = True
		print("Max. Q exited:   ", round(t,3), "s")
		print("Max. Q height:   ", round(alt/1000,3), "km")
		print("Max. Q speed:    ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("Max Q exit")

	# MECO
	if t>t_meco and not meco:
		meco = True
		print("MECO time:       ", round(t,3), "s")
		print("MECO height:     ", round(alt/1000,3), "km")
		print("MECO speed:      ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("MECO")

	# SEP 1
	if t>t_sep1 and not sep1:
		sep1 = True
		print("SEP1 time:       ", round(t,3), "s")
		print("SEP1 height:     ", round(alt/1000,3), "km")
		print("SEP1 speed:      ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("SEP1")

	# Second stage ignition
	if t>t_ssi and not ssi:
		ssi = True
		print("SSI time:        ", round(t,3), "s")
		print("SSI height:      ", round(alt/1000,3), "km")
		print("SSI speed:       ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("SSI")

	# SECO
	if t>t_seco and not seco:
		seco = True
		print("SECO time:       ", round(t,3), "s")
		print("SECO height:     ", round(alt/1000,3), "km")
		print("SECO speed:      ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("SECO")

	# SEP 2
	if t>t_sep2 and not sep2:
		sep2 = True
		print("SEP2 time:       ", round(t,3), "s")
		print("SEP2 height:     ", round(alt/1000,3), "km")
		print("SEP2 speed:      ", round(v,3), "m/s")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("SEP2")

	# Apogee
	if alt < alt_old and not apogee:
		apogee=True
		print("Apogee time:     ", round(t,3), "s")
		print("Apogee height:   ", round(alt/1000,3), "km")
		print("Apogee speed:    ", round(v,3), "m/s")
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event.append("Apogee")

	# Perigee
	if apogee and alt > alt_old and not perigee:
		perigee=True
		print("Perigee time:     ", round(t,3), "s")
		print("Perigee height:   ", round(alt/1000,3), "km")
		print("Perigee speed:    ", round(v,3), "m/s")
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event.append("Perigee")

	# Orbit
	if rp > R0+150000 and not orbit:
		orbit = True
		print("Orbit time:     ",round(t,3), "s")
		print("Apogee height:  ",round((ra-R0)/1000,3), "km")
		print("Perigee height: ",round((rp-R0)/1000,3), "km")
		list_event_x.append(x/1000)
		list_event_y.append(y/1000)
		list_event_a.append(alt/1000)
		list_event_d.append(d/1000)
		list_event.append("Orbit")

	# Once around
	if seco and theta_old<0 and theta>0:
		oncearound = True
		print("Once around time:     ",round(t,3), "s")

	# Impact
	if alt<=0 and liftoff:
		impact = True
		print("Impact time:     ",round(t,3), "s")
		print("Impact velocity: ",round(v,3), "m/s")
		print("Impact downrange distance: ",round(d/1000,3), "km")

	t = t + dt

print("Final FPA: %.2f deg" % (fpa/math.pi*180))
print("Final ha:  %.2f km" % ((ra-R0)/1000))
print("Final hp:  %.2f km" % ((rp-R0)/1000))
print("Final a:   %.2f km" % (a/1000))
print("Final e:   %.2f" % e)

# Planet
f, ax = plt.subplots(figsize=(6, 6))
circle = Circle((0,0),R0/1000,color='#000099')
ax.add_patch(circle)

# Orbit
plt.figure(1)
plt.plot(list_x,list_y,'b-')
plt.plot(list_event_x,list_event_y,'ro',markersize=5)
plt.xlabel('x')
plt.ylabel('y')
plt.grid(False)
plt.axis('scaled')
for i,txt in enumerate(list_event):
	ax.annotate(txt,(list_event_x[i],list_event_y[i]))

# Graphs
plt.figure(2)

plt.subplot(231)
plt.plot(list_t,list_alt,'b-',markersize=1)
plt.xlabel("Time / [s]")
plt.ylabel("Height / [km]")
plt.grid(True)

plt.subplot(232)
plt.plot(list_t,list_v,'b-',markersize=1)
plt.xlabel("Time / [s]")
plt.ylabel("Speed / [m/s]")
plt.grid(True)

plt.subplot(233)
plt.plot(list_t,list_F,'b-',markersize=1)
plt.xlabel("Time / [s]")
plt.ylabel("Thrust / [kN]")
plt.grid(True)

plt.subplot(234)
plt.plot(list_t,list_q,'b-',markersize=1)
plt.xlabel("Time / [s]")
plt.ylabel("Dynamic pressure / [Pa]")
plt.grid(True)

plt.subplot(235)
plt.plot(list_t,list_fpa,'b-',markersize=1)
plt.xlabel("Flight Path Angle / [deg]")
plt.ylabel("Time / [s]")
plt.grid(True)

plt.subplot(236)
plt.plot(list_v,list_alt,'b-',markersize=1)
plt.xlabel("Velocity / [m/s]")
plt.ylabel("Altitude / [km]")
plt.grid(True)

f, ax = plt.subplots()
ax.plot(list_d,list_alt)
plt.plot(list_event_d,list_event_a,'ro',markersize=5)
plt.xlabel("Downrange distance / [km]")
plt.ylabel("Altitude / [km]")
plt.axis('scaled')
plt.grid(True)
for i,txt in enumerate(list_event):
	ax.annotate(txt,(list_event_d[i],list_event_a[i]))

plt.show()
print("Done!")
