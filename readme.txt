This system simulates a system of spring forces using Forward Euler, Midpoint Method, Modified Midpoint Method, RK4 Method, Symplectic Euler, Backward Euler methods.


The file is saved in plain txt.

The movie is located at video.mp4.




My file format:  (.ptcls)
-- Control Separator --
(these are parameters, each line is one:)
useGravity
gravity
viscousDamping
springStiffness
springDamping
iterations
restitution
-- Below are Particles --
# of particles
index;pt.p0.x;pt.p0.y;pt.v0.x;pt.v0.y;pt.p.x;pt.p.y;pt.v.x;pt.v.y;pt.f.x;pt.f.y;pt.pinned
-- Below are Springs --
# of spring
p1.index; p2.index; restLength
