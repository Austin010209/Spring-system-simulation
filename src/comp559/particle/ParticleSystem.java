package comp559.particle;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.vecmath.Point2d;
import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;
import no.uib.cipr.matrix.sparse.CompDiagMatrix;
import mintools.parameters.BooleanParameter;
import mintools.parameters.DoubleParameter;
import mintools.parameters.IntParameter;
import mintools.swing.VerticalFlowPanel;
import mintools.viewer.SceneGraphNode;


/**
 * Implementation of a simple particle system
 * @author kry
 */
public class ParticleSystem implements SceneGraphNode, Function, Filter {
    
    private List<Particle> particles = new LinkedList<Particle>();
    
    private List<Spring> springs = new LinkedList<Spring>();
    
    /**
     * Creates an empty particle system
     */
    public ParticleSystem() {
        // do nothing
    }

    /**
     * Creates one of a number of simple test systems.
     * @param which
     */
    public void createSystem( int which ) {
        
        if ( which == 1) {        
            Point2d p = new Point2d( 100, 100 );
            Vector2d d = new Vector2d( 20, 0 );            
            Particle p1, p2, p3, p4;
            p1 = new Particle( p.x - d.y, p.y + d.x, 0, 0 );
            p1.index = particles.size();
            particles.add( p1 );
            p2 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
            p2.index = particles.size();
            particles.add( p2 );
            springs.add( new Spring ( p1, p2 ) );           
            p1.pinned = true;
            p2.pinned = true;            
            p.add( d );
            p.add( d );                    
            int N = 10;
            for (int i = 1; i < N; i++ ) {                
                //d.set( 20*Math.cos(i*Math.PI/N), 20*Math.sin(i*Math.PI/N) );                
                p3 = new Particle( p.x - d.y, p.y + d.x, 0, 0 );
                p3.index = particles.size();
                particles.add( p3 );
                p4 = new Particle( p.x + d.y, p.y - d.x, 0, 0 );
                p4.index = particles.size();
                particles.add( p4 );
                springs.add( new Spring ( p3, p1 ) );
                springs.add( new Spring ( p3, p2 ) );
                springs.add( new Spring ( p4, p1 ) );
                springs.add( new Spring ( p4, p2 ) );
                springs.add( new Spring ( p4, p3 ) );
                p1 = p3;
                p2 = p4;                
                p.add( d );
                p.add( d );            
            }
        } else if ( which == 2) {
            Particle p1 = new Particle( 320, 100, 0, 0 );
            p1.index = particles.size();
            particles.add( p1 );
            Particle p2 = new Particle( 320, 200, 0, 0 );
            p2.index = particles.size();
            particles.add( p2 );
            p1.pinned = true;
            springs.add( new Spring( p1, p2 ) );
        } else if ( which == 3 ) {
            int ypos = 100;
            Particle p0 = null;
            Particle p1, p2;
            p1 = new Particle( 320, ypos, 0, 0 );
            p1.index = particles.size();
            p1.pinned = true;            
            particles.add( p1 );
            int N = 10;
            for ( int i = 0; i < N; i++ ) {
                ypos += 20;
                p2 = new Particle( 320, ypos, 0, 0 );
                p2.index = particles.size();
                particles.add( p2 );
                springs.add( new Spring( p1, p2 ) );                
                // Hum.. this is not great in comparison to a proper bending energy...
                // use Maple to generate some code though, as it is painful to write by hand! :(
                if ( p0 != null ) springs.add( new Spring( p2, p0 ) );
                p0 = p1;
                
                p1 = p2;
            }
        }
    }
    
    /**
     * Gets the particles in the system
     * @return the particle set
     */
    public List<Particle> getParticles() {
        return particles;
    }
    
    /**
     * Gets the springs in the system
     * @return the spring list
     */
    public List<Spring> getSprings() {
    	return springs;
    }
    
    /**
     * Resets the positions of all particles
     */
    public void resetParticles() {
        for ( Particle p : particles ) {
            p.reset();
        }
        time = 0;
    }
    
    /**
     * Deletes all particles
     */
    public void clearParticles() {
        particles.clear();
        springs.clear();
    }
    
    /**
     * Gets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void getPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle p : particles ) {
            phaseSpaceState[count++] = p.p.x;
            phaseSpaceState[count++] = p.p.y;
            phaseSpaceState[count++] = p.v.x;
            phaseSpaceState[count++] = p.v.y;
        }
    }
    
    /**
     * Gets the dimension of the phase space state
     * (particles * 2 dimensions * 2 for velocity and position)
     * @return dimension
     */
    public int getPhaseSpaceDim() {        
        return particles.size() * 4;
    }
    
    /**
     * Sets the phase space state of the particle system
     * @param phaseSpaceState
     */
    public void setPhaseSpace( double[] phaseSpaceState ) {
        int count = 0;
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                count += 4;
            } else {
                p.p.x = phaseSpaceState[count++];
                p.p.y = phaseSpaceState[count++];
                p.v.x = phaseSpaceState[count++];
                p.v.y = phaseSpaceState[count++];
            }
        }
    }
    
    /**
     * Fixes positions and velocities after a step to deal with collisions 
     */
    public void postStepFix() {
        for ( Particle p : particles ) {
            if ( p.pinned ) {
                p.v.set(0,0);
            }
        }
        // do wall collisions
        double r = restitution.getValue();
        for ( Particle p : particles ) {            
            if ( p.p.x <= 0 ) {
                p.p.x = 0;
                if ( p.v.x < 0 ) p.v.x = - p.v.x * r;
                if ( p.f.x < 0 ) p.f.x = 0;                
            }
            if ( p.p.x >= width ) {
                p.p.x = width;
                if (p.v.x > 0 ) p.v.x = - p.v.x * r;
                if (p.f.x > 0 ) p.f.x = 0;
            } 
            
            if ( p.p.y >= height ) {
                p.p.y = height;
                if ( p.v.y > 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y > 0 ) p.f.y = 0;
            } 
            if ( p.p.y <= 0 ) {
                p.p.y = 0;
                if ( p.v.y < 0 ) p.v.y = - p.v.y * r;
                if ( p.f.y < 0 ) p.f.y = 0;
            }
        }
    }
    
    /** Elapsed simulation time */
    public double time = 0;

    /** The explicit integrator to use, if not performing backward Euler implicit integration */
    public Integrator integrator;
    
    public double[] state = new double[1];
    public double[] stateOut = new double[1];

    // these get created in init() and are probably useful for Backward Euler computations
    private ConjugateGradientMTJ CG;
    private DenseMatrix A;
    private DenseMatrix dfdx;
    private DenseMatrix dfdv;
    private DenseVector deltaxdot;
    private DenseVector b;
    private DenseVector f;
    private DenseVector xdot;
    
    /**
     * Initializes the system 
     * Allocates the arrays and vectors necessary for the solve of the full system
     */
    public void init() {
        int N = particles.size();
        // create matrix and vectors for solve
        CG = new ConjugateGradientMTJ(2*N);
        CG.setFilter(this);
        A = new DenseMatrix(2*N, 2*N);
        dfdx = new DenseMatrix(2*N, 2*N);
        dfdv = new DenseMatrix(2*N, 2*N);
        deltaxdot = new DenseVector(2*N);
        b = new DenseVector(2*N);
        f = new DenseVector(2*N);
        xdot = new DenseVector(2*N);
    }
    
    /**
     * Fills in the provided vector with the particle velocities.
     * @param xd
     */
    private void getVelocities(DenseVector xd) {
        for ( Particle p : particles ) {
            int j = p.index * 2;
            if( p.pinned ) {
                xd.set( j, 0 );
                xd.set( j+1, 0 );
            } else {
                xd.set( j, p.v.x );
                xd.set( j+1, p.v.y );
            }
        }       
    }

    /**
     * Sets the velocities of the particles given a vector
     * @param xd
     */
    private void setVelocities(DenseVector xd) {
        for ( Particle p : particles ) {
            int j = p.index * 2;
            if( p.pinned ) {
                p.v.set(0,0);
            } else {
                p.v.x = xd.get(j);
                p.v.y = xd.get(j+1);
            }
        }
    }
    
    /**
     *  Evaluates derivatives for ODE integration.
     * @param t time 
     * @param p phase space state
     * @param dydt to be filled with the derivative
     */
    @Override
    public void derivs(double t, double[] p, double[] dpdt) {
        // set particle positions to given values
        setPhaseSpace( p );
        
        // TODO: Objective 2, for explicit integrators, compute forces, and accelerations, and set dpdt
        //compute the accelerations and put them in the dpdt!
        
//        for(Particle pt: particles) {
//        	pt.clearForce();
//        	if(useGravity.getValue()) {
//        		double g = gravity.getValue();
//        		Vector2d grav = new Vector2d();
//        		grav.x = 0;
//        		grav.y = pt.mass * g;
//        		pt.addForce(grav);
//        	}
//        	double dampingcoef = -1*viscousDamping.getValue();
//        	Vector2d dampingforce = pt.v;
//        	dampingforce.scale(dampingcoef);
//        	pt.addForce(dampingforce);
//        }
        
        
        
        //Mine 2nd version
        for (Particle pt : particles) {
            pt.clearForce();
            if (useGravity.getValue()) {
                pt.f.y = gravity.getValue() * pt.mass;
            }
            pt.f.x -= viscousDamping.getValue() * pt.v.x;
            pt.f.y -= viscousDamping.getValue() * pt.v.y;
        }      
        
        
        for(Spring s: springs) {
        	s.apply();
        }
        
        
        int count = 0;
        for(Particle pt: particles) {
        	dpdt[count++] = pt.v.x;
        	dpdt[count++] = pt.v.y;
        	dpdt[count++] = pt.f.x / pt.mass;
        	dpdt[count++] = pt.f.y / pt.mass;
        	
        }
        
    }
    
    /** Time in seconds that was necessary to advance the system */
    public double computeTime;
    private DenseMatrix Mass;
    
    /**
     * Advances the state of the system
     * @param elapsed
     */
    public void advanceTime( double elapsed ) {
        Spring.k = springStiffness.getValue();
        Spring.c = springDamping.getValue();
            
        int n = getPhaseSpaceDim();
        
        long now = System.nanoTime();        
        
        if ( explicit.getValue() ) {
            if ( n != state.length ) {
                state = new double[n];
                stateOut = new double[n];
            }
            // TODO: See explicit stepping here
            
            getPhaseSpace(state);  
            integrator.step( state, n, time, elapsed, stateOut, this);                
            setPhaseSpace(stateOut);
            
        } else {        
            if ( f == null || f.size() != n ) {
                init();
            }
            
            // TODO: Objective 8, your backward Euler implementation will go here!
            // Note that the init() method called above creates a bunch of very 
            // useful MTJ working variables for you, and the ConjugateGradientMTJ object.
            // Go look at that code now!
            if ( n != state.length ) {
	            state = new double[n];
	            stateOut = new double[n];
            }
            
            //[M - h^2 dfdx - h dfdv] \delta \dot{x} = 
            // h f_0 + h^2 dfdx \dot{x_0}

            getPhaseSpace(state);
            if( !implicitNewTon.getValue() ) {
            	implicitEuler(state, n, elapsed, stateOut);
            }
            else {
                implicitEuler_Newton(linesearch.getValue(), state, n, elapsed, stateOut);
            }
            setPhaseSpace(stateOut);
        }
        time = time + elapsed;
        postStepFix();
        computeTime = (System.nanoTime() - now) / 1e9;
    }
    
    public void implicitEuler(double[] p, int n, double h, double[] pout) {
    	DenseVector x0 = new DenseVector(n/2);
    	DenseVector xout = new DenseVector(n/2);
    	DenseVector v0 = new DenseVector(n/2);
    	DenseVector vout = new DenseVector(n/2);
     	DenseVector temp = new DenseVector(n/2);
    	
    	for(int i = 0; i < n; i +=4) {
    		x0.set(i/2, p[i]);
    		x0.set(i/2+1, p[i+1]);
    		v0.set(i/2, p[i+2]);
    		v0.set(i/2+1, p[i+3]);
    	}
    	setup();
    	setupMass();
    	
        A.set( Mass.add(dfdx.scale(-h*h).add(dfdv.scale(-h))) );
        dfdx.mult(v0, temp);
        b.set(f.scale(h).add(temp.scale(h*h)));
        
        CG.solve(A, b, deltaxdot, iterations.getValue());
        
        vout.set(v0.add(deltaxdot));
        temp.set(vout);
        xout.set(x0.add(temp.scale(h)));
//        temp.scale(h);
//        xout.set(x0.add(temp.scale(h)));
        
        for(int i = 0; i < n; i+=4) {
        	pout[i] = xout.get(i/2);
        	pout[i+1] = xout.get(i/2+1);
        	pout[i+2] = vout.get(i/2);
        	pout[i+3] = vout.get(i/2+1);
        }
    }
    
    
    public void implicitEuler_Newton(boolean linesearch, double[] p0, int n, double h, double[] pout) {
    	//p should be the current phase space, and after
    	//this method, phase space will certainly change, so
    	//it is free to move, and no need to store initial p
    	//just store initial x0 and x0 dot, and keep them untouched
    	//because this method needs them.
    	
    	double tol = 0.0001;
    	double sigma = 0.0001;
    	
    	DenseVector v0 = new DenseVector(n/2);
    	DenseVector v1 = new DenseVector(n/2);
    	double[] newstate = new double[n];
    	DenseVector deltav = new DenseVector(n/2);
    	DenseVector Fv = new DenseVector(n/2);
    	DenseMatrix J = new DenseMatrix(n/2, n/2);
    	
    	for(int i = 0; i < n; i+=4) {
    		v0.set(i/2, p0[i+2]);
    		v0.set(i/2+1, p0[i+3]);
    		v1.set(i/2, p0[i+2]);
    		v1.set(i/2+1, p0[i+3]);
    	}
    	
    	getnewstate(p0, n, h, v1, newstate); //the first argument should always be p0: we need the x_0
    	setPhaseSpace(newstate);
    	setup();
    	setupMass();
    	
    	F_fun(v0, v1, h, Fv);  //similarly the first argument should always be v0
    	double lenFv = Fv.norm(Norm.Two);    	
    	
    	while( lenFv > tol ) {
        	F_deri_fun(J, h);
        	J.solve(Fv, deltav);
        	deltav.scale(-1);
        	if(!linesearch) {
        		v1.add(deltav);
        		
        		//for next while loop
        		getnewstate(p0, n, h, v1, newstate);
            	setPhaseSpace(newstate);
            	setup();
            	F_fun(v0, v1, h, Fv);
        	}
        	else {
        		double alpha = 1;
        		
        		v1.add(alpha, deltav);
        		getnewstate(p0, n, h, v1, newstate);
            	setPhaseSpace(newstate);
            	setup();
            	F_fun(v0, v1, h, Fv);
            	
            	while( (lenFv - Fv.norm(Norm.Two)) < sigma*alpha*lenFv ) {
            		alpha *= 0.5;
            		
            		//be ready for next loop
            		v1.add(-alpha, deltav);  //v1 was v + alpha delta v, and now it should be v + 0.5 alpha delta v
            		getnewstate(p0, n, h, v1, newstate);
                	setPhaseSpace(newstate);
                	setup();
                	F_fun(v0, v1, h, Fv);
            	}            	
        	}
        	lenFv = Fv.norm(Norm.Two);
    	}
    	
    	
    	for(int i = 0; i < n; i+=4) {
    		pout[i+2] = v1.get(i/2);  //assign v1
    		pout[i+3] = v1.get(i/2+1);
    		pout[i] = p0[i] + h*pout[i+2];
    		pout[i+1] = p0[i+1] + h*pout[i+3];
    	}
    	
    }
    
    
    public void setupMass() {
    	int n = 2*particles.size();
    	Mass = new DenseMatrix(n, n);
        for(int i = 0; i < n; i++) {
        	Mass.set(i, i, particles.get(i/2).mass);
        }
    }
    
    
    //input of f
    public void getnewstate(double[] p, int n, double h, DenseVector xdot, double[] pout) {
    	for(int i = 0; i < n; i+=4) {
    		pout[i+2] = xdot.get(i/2);
    		pout[i+3] = xdot.get(i/2+1);
    		pout[i] = p[i] + h*pout[i+2];
    		pout[i+1] = p[i+1] + h*pout[i+3];
    	}
    }
    
    public void F_fun(DenseVector v0, DenseVector v1, double h, DenseVector ans) {
    	Mass.solve(f, ans);
    	ans.scale(h);
    	ans.add(v0);
    	ans.add(-1, v1);
    }
    
    public void F_deri_fun(DenseMatrix ans, double h){
    	int n = particles.size()*2;
    	
    	//a potential solve to the problem
    	DenseMatrix temp = new DenseMatrix(n,n);
    	DenseMatrix temp1 = new DenseMatrix(n,n);
    	temp.set(dfdv);
    	
    	//
//    	for(int i = 0; i < n; i+=2) {
//    		temp.add(i, i, -viscousDamping.getValue());
//    		temp.add(i+1, i+1, -viscousDamping.getValue());
//    		temp.add(i+1, i+1, -gravity.getValue()); 
//    	}
    	//
    	
    	temp1.add(dfdx);
    	temp1.scale(h);
    	temp.add(temp1);
    	Mass.solve(temp, ans);
    	ans.scale(h);
    	
    	
    	DenseMatrix I = new DenseMatrix(n, n);
    	for(int i = 0; i < n; i++) {
    		I.set(i, i, -1);
    	}
    	ans.add(I);
    }
    
    public void setup() {
    	dfdx.zero();
        dfdv.zero();
        f.zero();
    	for(Spring s: springs) {
        	s.addForce(f);
        	s.addDfdx(dfdx);
        	s.addDfdv(dfdv);
        }
        
        for (Particle pt : particles)
        {
            f.add(pt.index * 2, -viscousDamping.getValue() * pt.v.x);
            f.add(pt.index * 2 + 1, -viscousDamping.getValue() * pt.v.y);
            f.add(pt.index * 2 + 1, -gravity.getValue() * pt.v.y);
        }        
    }
    
    
    
    @Override
    public void filter(Vector v) {
        for ( Particle p : particles ) {
            if ( !p.pinned ) continue;
            v.set( p.index*2+0, 0 );
            v.set( p.index*2+1, 0 );
        }
    }

    /**
     * Creates a new particle and adds it to the system
     * @param x
     * @param y
     * @param vx
     * @param vy
     * @return the new particle
     */
    public Particle createParticle( double x, double y, double vx, double vy ) {
        Particle p = new Particle( x, y, vx, vy );
        p.index = particles.size();
        particles.add( p );
        return p;
    }
    
    public void remove( Particle p ) {
    	for ( Spring s : p.springs ) {
    		Particle other = s.p1 == p ? s.p2 : s.p1; 
    		other.springs.remove( s );
    		springs.remove( s );
    	}
    	p.springs.clear(); // not really necessary
    	particles.remove( p );
    	// reset indices of each particle :(
    	for ( int i = 0 ; i < particles.size(); i++ ) {
    		particles.get(i).index = i;
    	}
    }
    
    /**
     * Creates a new spring between two particles and adds it to the system.
     * @param p1
     * @param p2
     * @return the new spring
     */
    public Spring createSpring( Particle p1, Particle p2 ) {
        Spring s = new Spring( p1, p2 ); 
        springs.add( s );         
        return s;
    }
    
    /**
     * Removes a spring between p1 and p2 if it exists, does nothing otherwise
     * @param p1
     * @param p2
     * @return true if the spring was found and removed
     */
    public boolean removeSpring( Particle p1, Particle p2 ) {
    	Spring found = null;
    	for ( Spring s : springs ) {
    		if ( ( s.p1 == p1 && s.p2 == p2 ) || ( s.p1 == p2 && s.p2 == p1 ) ) {
    			found = s;
    			break;
    		}
    	}
    	if ( found != null ) {
    		found.p1.springs.remove(found);
    		found.p2.springs.remove(found);
    		springs.remove(found);
			return true;
    	}
    	return false;
    }
    
    @Override
    public void init(GLAutoDrawable drawable) {
        // do nothing
    }

    private int height;
    private int width;

    @Override
    public void display(GLAutoDrawable drawable) {
        GL2 gl = drawable.getGL().getGL2();

        // update the width and the height for wall collisions
        height = drawable.getSurfaceHeight();
        width = drawable.getSurfaceWidth();
        
        gl.glPointSize( 10 );
        gl.glBegin( GL.GL_POINTS );
        for ( Particle p : particles ) {
            double alpha = 0.5;
            if ( p.pinned ) {
                gl.glColor4d( 1, 0, 0, alpha );
            } else {
                gl.glColor4d( p.color.x, p.color.y, p.color.z, alpha );
            }
            gl.glVertex2d( p.p.x, p.p.y );
        }
        gl.glEnd();
        
        gl.glColor4d(0,.5,.5,.5);
        gl.glLineWidth(2f);
        gl.glBegin( GL.GL_LINES );
        for (Spring s : springs) {
            gl.glVertex2d( s.p1.p.x, s.p1.p.y );
            gl.glVertex2d( s.p2.p.x, s.p2.p.y );
        }
        gl.glEnd();
    }
    
    public BooleanParameter useGravity = new BooleanParameter( "use gravity", true );
    public DoubleParameter gravity = new DoubleParameter( "gravity", 9.8, 0.01, 1000 );
    public DoubleParameter springStiffness = new DoubleParameter( "spring stiffness", 100, 0, 10000 );
    public DoubleParameter springDamping = new DoubleParameter( "spring damping", 0, 0, 50 );
    public DoubleParameter viscousDamping = new DoubleParameter( "viscous damping", 0, 0, 10 );
    public DoubleParameter restitution = new DoubleParameter( "r", 0, 0, 1 );
    public JTextArea comments = new JTextArea("enter comments in control panel");
    public IntParameter iterations = new IntParameter( "iterations", 100, 1, 100 );
    /** controls weather explicit or implicit integration is used */
    public BooleanParameter explicit = new BooleanParameter( "explicit", true );
    //
    public BooleanParameter implicitNewTon = new BooleanParameter("newton", false);
    public BooleanParameter linesearch = new BooleanParameter("linesearch", false);

    
    @Override
    public JPanel getControls() {
        VerticalFlowPanel vfp = new VerticalFlowPanel();
        vfp.add( comments );
        vfp.add( useGravity.getControls() );
        vfp.add( gravity.getSliderControls(true) );
        vfp.add( springStiffness.getSliderControls(false) );
        vfp.add( springDamping.getSliderControls(false) );
        vfp.add( viscousDamping.getSliderControls(false) );
        vfp.add( restitution.getSliderControls(false) );
        vfp.add( iterations.getSliderControls() );
        vfp.add( explicit.getControls() );
        vfp.add( implicitNewTon.getControls() );
        vfp.add( linesearch.getControls() );
        return vfp.getPanel();        
    }
    
    @Override
    public String toString() {
        // TODO: Add your name below
        String ret = "Austin Zhang\n" +
                     comments.getText() + "\n" +
                     "particles = " + particles.size() + "\n";
        if ( explicit.getValue() ) {
            ret += "integrator = " + integrator.getName() + "\n";
        } else {
            ret += "integrator = Backward Euler\n";
        }
        ret += "k = " + springStiffness.getValue() + "\n" +
               "c = " + springDamping.getValue() + "\n" +
               "b = " + viscousDamping.getValue() +"\n" + 
               "time = " + time;
        return ret;
    }
    
}
