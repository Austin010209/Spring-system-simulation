package comp559.particle;

import javax.vecmath.Vector2d;

import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.Matrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.Vector.Norm;

/**
 * Spring class for 599 assignment 1
 * @author kry
 */
public class Spring {

    Particle p1 = null;
    Particle p2 = null;
    
    /** Spring stiffness, sometimes written k_s in equations */
    public static double k = 1;
    /** Spring damping (along spring direction), sometimes written k_d in equations */
    public static double c = 1;
    /** Rest length of this spring */
    double l0 = 0;
    
    private boolean Mouse_is_Spring = false;
    private double M_K = 1;
    private double M_C = 1;
    
    /**
     * Creates a spring between two particles
     * @param p1
     * @param p2
     */
    public Spring( Particle pp1, Particle pp2 ) {
        p1 = pp1;
        p2 = pp2;
        recomputeRestLength();
        p1.springs.add(this);
        p2.springs.add(this);
    }
    
    /**
     * Computes and sets the rest length based on the original position of the two particles 
     */
    public void recomputeRestLength() {
        l0 = p1.p0.distance( p2.p0 );
    }
    
    
    public void setMouseSpring(double pM_K, double pM_C) {
        M_C = pM_C;
        M_K = pM_K;
        Mouse_is_Spring = true;
    }
    
    /**
     * Applies the spring force by adding a force to each particle
     */
    public void apply() {
        // TODO: Objective 1, FINISH THIS CODE!
    	Vector2d dif = new Vector2d();
    	dif.sub(p1.p, p2.p);  //A-B
    	double len = dif.length();
    	dif.normalize();  //direction

    	dif.scale(-k*(len-l0));
        p1.addForce(dif);
        dif.scale(-1);
        p2.addForce(dif);
        dif.scale(-1);
        
        dif = new Vector2d();
    	dif.sub(p1.p, p2.p);
        dif.normalize();
        Vector2d damp = new Vector2d(dif);
        
        Vector2d v = new Vector2d();
        v.sub(p2.v, p1.v);
        damp.scale(-c * damp.dot(v)); 
        
        p1.addForce(damp);
        damp.scale(-1);
        p2.addForce(damp);
        
    }
   
    /** TODO: the functions below are for the backwards Euler solver */
    
    /**
     * Computes the force and adds it to the appropriate components of the force vector.
     * (This function is something you might use for a backward Euler integrator)
     * @param f
     */
    public void addForce( Vector f ) {
        // TODO: Objective 8, FINISH THIS CODE for backward Euler method (probably very simlar to what you did above)
    	Vector2d springforce = new Vector2d();
    	springforce.sub(p1.p, p2.p);
    	double len = springforce.length();
    	springforce.normalize();
    	
    	Vector2d dampingforce = new Vector2d(springforce);
    	
    	springforce.scale(-k*(len - l0));
    	
    	Vector2d velocity = new Vector2d();
    	velocity.sub(p1.v, p2.v);
    	dampingforce.scale( -c * dampingforce.dot(velocity) );
    	
    	Vector2d force = new Vector2d();
    	force.add(springforce, dampingforce);
    	
    	f.add(p1.index*2, force.x);
    	f.add(p1.index*2+1, force.y);
    	
    	force.scale(-1);
    	f.add(p2.index*2, force.x);
    	f.add(p2.index*2+1, force.y);
    }
    
    /**
     * Adds this springs contribution to the stiffness matrix
     * @param dfdx
     */
    public void addDfdx( Matrix dfdx ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward euler integration
        DenseVector pos = new DenseVector(2);
        pos.set(0, p1.p.x-p2.p.x);
        pos.set(1, p1.p.y-p2.p.y);

        
//        double theta = Math.atan(pos.y/pos.x);
//        if(Math.signum(pos.y/pos.x)==-1)theta += Math.PI;
//        double cos = Math.cos(theta);
//        double sin = Math.sin(theta);
//        double[][] mat = {{cos, -sin}, {sin, cos}};
//        Matrix R = new DenseMatrix(mat);
//        double[][] mat1 = {{cos, sin}, {-sin, cos}};
//        Matrix RT = new DenseMatrix(mat1);
//        double[][] mat2 = {{k, 0}, {0, 0}};
//        Matrix K = new DenseMatrix(mat2);
//        Matrix dFada = new DenseMatrix(mat2);
//        R.mult(-1, K, dFada);
//        dFada.mult(RT, dFada);
        
        // -k/ ||l||^2 * l*lT
        double len = pos.norm(Norm.Two);
        Matrix dFada = new DenseMatrix(2,2);
        dFada.rank1(pos);
        dFada.scale(-k/(len*len));
        
        
        
        int[][] order = {{0,0}, {0,1}, {1,0}, {1,1}};
        int posA = p1.index*2;
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdx.add(posA+a, posA+b, dFada.get(a, b));
        }
        int posB = p2.index*2;
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdx.add(posB+a, posB+b, dFada.get(a, b));
        }
        dFada.scale(-1);
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdx.add(posA+a, posB+b, dFada.get(a, b));
        }
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdx.add(posB+a, posA+b, dFada.get(a, b));
        }
    }   
 
    /**
     * Adds this springs damping contribution to the implicit damping matrix
     * @param dfdv
     */
    public void addDfdv( Matrix dfdv ) {
        // TODO: Objective 8, FINISH THIS CODE... necessary for backward Euler integration
        // The thing only depends on velocity: that is viscous damping!(not this one)
    	
    	Vector2d direction = new Vector2d();
    	direction.sub(p1.p, p2.p);
    	direction.normalize();
    	
    	DenseVector posdif = new DenseVector(2);
    	DenseMatrix dFdAv = new DenseMatrix(2,2);
    	posdif.set(0, direction.x);
    	posdif.set(1, direction.y);
    	dFdAv.rank1(posdif);
    	dFdAv.scale(-c);
    	
    	
    	int[][] order = {{0,0}, {0,1}, {1,0}, {1,1}};
        int posA = p1.index*2;
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdv.add(posA+a, posA+b, dFdAv.get(a, b));
        }
        int posB = p2.index*2;
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdv.add(posB+a, posB+b, dFdAv.get(a, b));
        }
        dFdAv.scale(-1);
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdv.add(posA+a, posB+b, dFdAv.get(a, b));
        }
        for(int j = 0; j < order.length; j++) {
        	int a = order[j][0];
        	int b = order[j][1];
        	dfdv.add(posB+a, posA+b, dFdAv.get(a, b));
        }
    	
    } 
    
}
