package comp559.particle;

public class ModifiedMidpoint implements Integrator {

	private double[] tmp;
	
    @Override
    public String getName() {
        return "modified midpoint";
    }

    @Override
    public void step(double[] p, int n, double t, double h, double[] pout, Function derivs) {
    	// TODO: Objective 5, implmement the modified midpoint (2/3) method.
    	// see also efficient memory management suggestion in provided code for the Midpoint method.
    	
    	if ( tmp == null || tmp.length != n ) {
            tmp = new double[n];
    	}
        derivs.derivs(t, p, tmp);  //dpdt
        
    	for(int i = 0; i < n; i++) {
    		tmp[i] = p[i] + 2.0/3*h*tmp[i];   //half_p   dpdt
    	}
        derivs.derivs(t + 2.0/3.0*h, tmp, tmp);  //half_p, half_dpdt
        
        for(int i = 0; i < n; i++) {
        	pout[i] = p[i] + h*tmp[i];  //half_dpdt
        }
    }

}
