package comp559.particle;

public class RK4 implements Integrator {
    
	private double[] tmp;
	private double[] increment;
	
    @Override
    public String getName() {
        return "RK4";
    }

    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
        // TODO: Objective 6, implement the RK4 integration method
    	// see also efficient memory management suggestion in provided code for the Midpoint method.
    	
    	if ( tmp == null || tmp.length != n ) {
            tmp = new double[n];
    	}
    	if ( increment == null || increment.length != n ) {
    		increment = new double[n];
    	}
    	
    	derivs.derivs(t, p, tmp);  //k1
    	increment = tmp;
    	for(int i = 0; i < n; i++) {
    		tmp[i] = p[i] + 0.5*h * tmp[i];
    	}
    	derivs.derivs(t, tmp, tmp);  //k2
    	for(int i = 0; i < n; i++) {
    		increment[i] += 2*tmp[i];
    	}
    	
    	for(int i = 0; i < n; i++) {
    		tmp[i] = p[i] + 0.5*h * tmp[i];
    	}
    	derivs.derivs(t, tmp, tmp);  //k3
    	for(int i = 0; i < n; i++) {
    		increment[i] += 2*tmp[i];
    	}
    	
    	for(int i = 0; i < n; i++) {
    		tmp[i] = p[i] + h * tmp[i];
    	}
    	derivs.derivs(t, tmp, tmp);  //k4
    	for(int i = 0; i < n; i++) {
    		increment[i] += tmp[i];
    	}
    	
    	for(int i = 0; i < n; i++) {
    		pout[i] = p[i] + h/6.0*increment[i];
    	}
    }
}
