package optimists.hbv;

public class HBV
{

	// --------------------------------------------------------------------------------------------
	// Main
	// --------------------------------------------------------------------------------------------
	
	public static void main(String[] args)
	{
		double tt				= 0.0;
		double dd				= 0.29;
		double fc				= 165.0;
		double beta				= 1.8;
		double pelt				= 7.42;
		double tlt				= -3.14;
		double c				= 0.078;
		double wp				= 22.5;
		double kPerc			= 0.019;
		double k1				= 0.013;
		double k2				= 0.023;
		
		double s_init			= 0.0;
		double sm1_init			= 10.0;
		double sm2_init			= 10.0;
		
		double areaMult			= 3031/(24*3.6);
		
		double[] precip			= {0, 0.11, 1.55, 0, 0, 0, 0, 0, 0, 0.27, 9.17, 0, 0, 0, 0, 0,
				13.55, 0, 0, 0, 0, 2.31, 1.6, 2.63, 0.6, 6.15, 0, 0, 0, 0.11, 0.13, 0, 0, 0, 0, 0,
				0, 0, 0.11, 9.83, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.31, 13.09, 0.04, 41.87, 0.17, 0, 0,
				0, 0, 5.97, 35.56, 5.73, 0, 0, 0, 0, 0, 29.26, 7.49, 0, 0, 1.43, 2.61, 0, 0, 0, 0,
				0, 0, 0, 0, 0, 0, 10.35, 0, 8.29, 0, 0, 0, 0, 0, 0, 0, 6.44, 0, 1.11, 0, 0, 0, 0, 
				0, 13.46, 0.02, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0, 0, 5.67, 4.88, 8.84, 2.53,
				34.93, 24.39, 0, 33.95, 7.87, 0, 0, 0, 1.32, 0, 9.13, 0, 0.41, 1.21, 6.78, 0, 0.12,
				0, 0.28, 0, 0, 0, 0, 0, 0, 7.42, 13.27, 53.62, 0.96, 0.01, 22.47, 0, 0, 0, 0, 0, 0,
				8.63, 6.25, 0, 0.42, 13.43, 8.29, 0.79, 0, 0, 0, 2.42, 0.02, 0, 1.93, 0.83, 0.17, 
				1.49, 1.09, 0.93, 0, 0.23, 0, 0, 0, 0, 3.5, 0.13};
		
		double[] temp			= {-1.51, -0.58, -3.44, -3.15, -1.76, 1.48, -8.76, -8.87, -4.4, 
				3.74, 6.34, 2.31, -0.87, -2.67, -7.09, -4.55, -6.84, -10.75, -10.49, -4, -3.69, 
				-0.71, 5.01, 7.73, 11.02, 11.89, -0.44, -3.76, -5.51, -4.73, -7.88, -10.19, -8.82,
				-5.01, 2.58, 1.37, -2.28, 0.11, 0.08, -11.67, -11.14, -8.21, -2.32, -6.11, -4.52,
				1.27, -0.27, 1.82, 7.05, 10.98, 3.84, 3.21, 4.53, -2.16, -5.76, -2.93, -8.81,
				-7.05, 0.19, 3.61, -0.87, 0.42, 3.96, 6.71, 9.76, 12.84, 2.85, -1.54, -3.44, 1.49,
				4.64, 6.12, 2.82, 9.55, 2.99, 8.76, 11.36, 11.66, 15.81, 7.97, 2.58, 13.16, 13.02,
				3.71, 9.43, 7.6, 0.85, 3.03, 2.9, 1.41, 7.05, 11.18, 5.75, 5.72, 6.46, -0.56,
				0.16, 7.44, 11.07, 16.9, 15.69, 5.27, 5.8, 12.42, 16, 4.06, 5.5, 10.08, 13.59,
				13.73, 13.97, 14.01, 14.54, 15.27, 18.65, 18.98, 16.13, 7.87, 7.78, 6.46, 8.01, 
				8.83, 9.18, 10.84, 13.91, 15.58, 16.85, 9.04, 13.68, 15.33, 16.54, 18.93, 18.01, 
				16.99, 16.59, 14.65, 15.97, 15.7, 12.79, 12.89, 12.98, 14.7, 15.99, 17.97, 17.59,
				17.68, 15.69, 16.08, 16.58, 17.06, 17.66, 18.92, 19.21, 20.11, 20.54, 21.62, 18.83,
				19.63, 24.84, 24.54, 20.31, 19.25, 18.68, 22.08, 23.33, 22.59, 20.79, 21.75, 20.34,
				21.05, 21.67, 20.9, 20.66, 22.33, 18.66, 21.91, 23.72, 24.98, 22.74, 24.27, 22.64,
				23.68};
		
		HBVElement watershed	= new HBVElement(tt, dd, fc, beta, pelt, tlt, c, wp, kPerc, k1,
													k2, s_init, sm1_init, sm2_init);
		
		double[] q				= runExtendedSimulation(watershed, precip, temp);
		
		System.out.println("Streamflow (m3/s)");
		for (int t = 0; t < q.length; t++)
			System.out.println(q[t]*areaMult);
	}
	
	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Error message: inputs arrays do not have equal sizes
	 */
	public final static String ERR_INPUT_SIZE = "Sizes of input arrays do not match: %1$f, %2$f";
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Runs a extended time simulation using a lumped HBV model
	 * @param watershed		Contains the parameters and initial conditions of the model
	 * @param precipitation	The time series of precipitation depth for each time step
	 * @param temperature	The time series of temperature for each time step
	 * @return				The time series of output water depth for each time step 
	 */
	public static double[] runExtendedSimulation(HBVElement watershed, double[] precipitation,
													double[] temperature)
	{
		int nPrec		= precipitation.length;
		int nTemp		= temperature.length;
		if (nPrec != nTemp)
			throw new IllegalArgumentException(String.format(ERR_INPUT_SIZE, nPrec, nTemp));
		
		double[] q		= new double[nPrec];
		for (int t = 0; t < nPrec; t++)
		{
			double prec	= precipitation[t];
			double temp	= temperature[t];
			watershed.simulateTimeStep(prec, temp);
			q[t]		= watershed.getQ();
		}
		return q;
	}

}
