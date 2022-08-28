package comp559.particle;

public class Midpoint implements Integrator {

    @Override
    public String getName() {
        return "midpoint";
    }

    private double[] tmp;
    
    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 4, implement midpoint method

    	// You will probably want a temporary array in this method and perhaps
    	// multiple temporary arrays in other higher order explicit integrators.
    	// Avoid thrashing memory by reallocating only when necessary.
    	if ( tmp == null || tmp.length != n ) {
            tmp = new double[n];
    	}
        derivs.derivs(t, p, tmp);  //dpdt
        
    	for(int i = 0; i < n; i++) {
    		tmp[i] = p[i] + 0.5*h*tmp[i];   //half_p   dpdt
    	}
        derivs.derivs(t + 0.5*h, tmp, tmp);  //half_p, half_dpdt
        
        for(int i = 0; i < n; i++) {
        	pout[i] = p[i] + h*tmp[i];  //half_dpdt
        }
    }

}
