package optimists.hbv;

/**
 * Represents a watershed element of the HBV (Hydrologiska Byråns Vattenbalansavdelning)
 * hydrologic modeling engine. Allows to perform step-wise hydrologic calculations (including soil
 * moisture, evapotranspiration, and streamflow) provided precipitation and air temperature inputs.
 * Any units for depths [L], temperatures [H], and the time step [T] can be used, but they need to
 * be consistent throughout the parameters and forcings. 
 * Based on:  A. Aghakouchak and E. Habib, “Application of a conceptual hydrologic model in
 * teaching hydrologic processes,” Int. J. Eng. Educ., vol. 26, no. 4, pp. 963–973, 2010.
 * @author Felipe Hernández
 */
public class HBVElement
{

	// --------------------------------------------------------------------------------------------
	// Parameters
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Snowmelt threshold temperature: snow melts if the temperature is higher than this value;
	 * precipitation accumulates as snow otherwise [H]
	 */
	private double tt;
	
	/**
	 * Snowmelt "degree-day" factor: the rate at which {@link #s} decreases per degree of 
	 * temperature above {@link #tt} and per time step [LH<sup>-1</sup>T<sup>-1</sup>]
	 */
	private double dd;
	
	/**
	 * Field capacity: the maximum value for {@link #sm1} [L]
	 */
	private double fc;
	
	/**
	 * Shape coefficient: exponent of the effective precipitation generated based on the ratio
	 * between {@link #sm1} and {@link #fc} [-]
	 */
	private double beta;
	
	/**
	 * Long-term mean potential rate of evapotranspiration per time step for the current season.
	 * Used to estimate the current potential evapotranspiration depending on the temperature.
	 * [LT<sup>-1</sup>]
	 */
	private double pelt;
	
	/**
	 * Long-term mean temperature in the day for the current season. Used to estimate the 
	 * potential evapotranspiration, which is based on the difference between this value and the
	 * current temperature. [H]
	 */
	private double tlt;
	
	/**
	 * Potential evapotranspiration coefficient: to compute the potential rate of
	 * evapotranspiration depending on the difference between the temperature and {@link #tlt}
	 * [H<sup>-1</sup>]
	 */
	private double c;
	
	/**
	 * Wilting point: value of {@link #sm1} below which {@link #e} becomes limited due to soil
	 * dryness [L]
	 */
	private double wp;
	
	/**
	 * Percolation coefficient: coefficient to compute {@link #qPerc} based on {@link #sm1}
	 * [T<sup>-1</sup>]
	 */
	private double kPerc;
	
	/**
	 * Near surface runoff storage coefficient: coefficient to compute {@link #q1} based on
	 * {@link #sm1} [T<sup>-1</sup>]
	 */
	private double k1;
	
	/**
	 * Baseflow storage coefficient: coefficient to compute {@link #q2} based on {@link #sm2}
	 * [T<sup>-1</sup>]
	 */
	private double k2;
	
	// --------------------------------------------------------------------------------------------
	// State variables
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Depth of the water-equivalent layer of the snow pack [L]
	 */
	private double s;
	
	/**
	 * Depth of the water contained as soil moisture in the top soil layer. Soil layers do not have
	 * a specific depth. [L]
	 */
	private double sm1;
	
	/**
	 * Depth of the water contained as soil moisture in the bottom soil layer. Soil layers do not
	 * have a specific depth. [L]
	 */
	private double sm2;
	
	// --------------------------------------------------------------------------------------------
	// Output
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Snowmelt: equivalent water depth that melted from {@link #s} [L]
	 */
	private double sm;
	
	/**
	 * Effective precipitation: amount of precipitation that directly became runoff [L]
	 */
	private double pEff;
	
	/**
	 * Actual potential evapotranspiration: depth of equivalent water that can be released to the
	 * atmosphere in ideal conditions [L]
	 */
	private double pe;
	
	/**
	 * Evapotranspiration: depth of equivalent water that was released from {@link #sm1} to the
	 * atmosphere [L]
	 */
	private double e;
	
	/**
	 * Percolation: amount of {@link #sm1} that was transferred to {@link #sm2} [L]
	 */
	private double qPerc;
	
	/**
	 * Near-surface runoff: amount of {@link #sm1} released to the drainage network [L]
	 */
	private double q1;
	
	/**
	 * Baseflow: amount of {@link #sm2} released to the drainage network [L]
	 */
	private double q2;
	
	/**
	 * Streamflow: total water depth that was released to the drainage network [L]
	 */
	private double q;

	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param tt		{@link #tt}
	 * @param dd		{@link #dd}
	 * @param fc		{@link #fc}
	 * @param beta		{@link #beta}
	 * @param pelt		{@link #pelt}
	 * @param tlt		{@link #tlt}
	 * @param c			{@link #c}
	 * @param wp		{@link #wp}
	 * @param kPerc		{@link #kPerc}
	 * @param k1		{@link #k1}
	 * @param k2		{@link #k2}
	 * @param s_init	initial value for {@link #s}
	 * @param sm1_init	initial value for {@link #sm1}
	 * @param sm2_init	initial value for {@link #sm2}
	 */
	public HBVElement(double tt, double dd, double fc, double beta, double pelt, double tlt,
			double c, double wp, double kPerc, double k1, double k2, double s_init,
			double sm1_init, double sm2_init)
	{
		this.tt		= tt;
		this.dd		= dd;
		this.fc		= fc;
		this.beta	= beta;
		this.pelt	= pelt;
		this.tlt	= tlt;
		this.c		= c;
		this.wp		= wp;
		this.kPerc	= kPerc;
		this.k1		= k1;
		this.k2		= k2;
		
		this.s		= Math.max(0, s_init	);
		this.sm1	= Math.max(0, sm1_init	);
		this.sm2	= Math.max(0, sm2_init	);
		
		this.pEff	= Double.NaN;
		this.pe		= Double.NaN;
		this.e		= Double.NaN;
		this.qPerc	= Double.NaN;
		this.q1		= Double.NaN;
		this.q2		= Double.NaN;
		this.q		= Double.NaN;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * @return {@link #tt}
	 */
	public double getTt()
	{
		return tt;
	}

	/**
	 * @param tt {@link #tt}
	 */
	public void setTt(double tt)
	{
		this.tt = tt;
	}

	/**
	 * @return {@link #dd}
	 */
	public double getDd()
	{
		return dd;
	}

	/**
	 * @param dd {@link #dd}
	 */
	public void setDd(double dd)
	{
		this.dd = dd;
	}

	/**
	 * @return {@link #fc}
	 */
	public double getFc()
	{
		return fc;
	}

	/**
	 * @param fc {@link #fc}
	 */
	public void setFc(double fc)
	{
		this.fc = fc;
	}

	/**
	 * @return {@link #beta}
	 */
	public double getBeta()
	{
		return beta;
	}

	/**
	 * @param beta {@link #beta}
	 */
	public void setBeta(double beta)
	{
		this.beta = beta;
	}

	/**
	 * @return {@link #pelt}
	 */
	public double getPelt()
	{
		return pelt;
	}

	/**
	 * @param pelt {@link #pelt}
	 */
	public void setPelt(double pelt)
	{
		this.pelt = pelt;
	}

	/**
	 * @return {@link #tlt}
	 */
	public double getTlt()
	{
		return tlt;
	}

	/**
	 * @param tlt {@link #tlt}
	 */
	public void setTlt(double tlt)
	{
		this.tlt = tlt;
	}

	/**
	 * @return {@link #c}
	 */
	public double getC()
	{
		return c;
	}

	/**
	 * @param c {@link #c}
	 */
	public void setC(double c)
	{
		this.c = c;
	}

	/**
	 * @return {@link #wp}
	 */
	public double getWp()
	{
		return wp;
	}

	/**
	 * @param wp {@link #wp}
	 */
	public void setWp(double wp)
	{
		this.wp = wp;
	}

	/**
	 * @return {@link #kPerc}
	 */
	public double getkPerc()
	{
		return kPerc;
	}

	/**
	 * @param kPerc {@link #kPerc}
	 */
	public void setkPerc(double kPerc)
	{
		this.kPerc = kPerc;
	}

	/**
	 * @return {@link #k1}
	 */
	public double getK1()
	{
		return k1;
	}

	/**
	 * @param k1 {@link #k1}
	 */
	public void setK1(double k1)
	{
		this.k1 = k1;
	}

	/**
	 * @return {@link #k2}
	 */
	public double getK2()
	{
		return k2;
	}

	/**
	 * @param k2 {@link #k2}
	 */
	public void setK2(double k2)
	{
		this.k2 = k2;
	}

	/**
	 * @return {@link #s}
	 */
	public double getS()
	{
		return s;
	}

	/**
	 * @param s {@link #s}
	 */
	public void setS(double s)
	{
		this.s = s;
	}

	/**
	 * @return {@link #sm1}
	 */
	public double getSm1()
	{
		return sm1;
	}

	/**
	 * @param sm1 {@link #sm1}
	 */
	public void setSm1(double sm1)
	{
		this.sm1 = sm1;
	}

	/**
	 * @return {@link #sm2}
	 */
	public double getSm2()
	{
		return sm2;
	}

	/**
	 * @param sm2 {@link #sm2}
	 */
	public void setSm2(double sm2)
	{
		this.sm2 = sm2;
	}

	/**
	 * @return {@link #sm}
	 */
	public double getSm()
	{
		return sm;
	}

	/**
	 * @return {@link #pEff}
	 */
	public double getpEff()
	{
		return pEff;
	}

	/**
	 * @return {@link #pe}
	 */
	public double getPe()
	{
		return pe;
	}

	/**
	 * @return {@link #e}
	 */
	public double getE()
	{
		return e;
	}

	/**
	 * @return {@link #qPerc}
	 */
	public double getqPerc()
	{
		return qPerc;
	}

	/**
	 * @return {@link #q1}
	 */
	public double getQ1()
	{
		return q1;
	}

	/**
	 * @return {@link #q2}
	 */
	public double getQ2()
	{
		return q2;
	}

	/**
	 * @return {@link #q}
	 */
	public double getQ()
	{
		return q;
	}
	
	/**
	 * Advances the hydrologic simulation one time step. This computes the values for the updated
	 * state variables and the output variables.
	 * @param p	Precipitation: the depth of liquid and solid water that precipitates from the
	 * 			atmosphere during the time step [L]
	 * @param t	Air temperature during the time step [H]
	 */
	public void simulateTimeStep(double p, double t)
	{
		double p_liquid	= t < tt ? 0.0 : p;
		double p_solid	= p - p_liquid;
		
		sm				= s > 0.0 ? (t > tt ? Math.max(s, dd*(t - tt)) : 0.0) : 0.0;
		s				= Math.max(0.0, s + p_solid - sm);
		pEff			= Math.pow(sm1/fc, beta)*(p_liquid + sm);
		pe				= (1.0 + c*(t - tlt))*pelt;
		qPerc			= sm1*Math.min(1.0, kPerc);
		q1				= Math.min(sm1 - qPerc, sm1*k1);
		q2				= Math.min(sm2, sm2*k2);
		e				= Math.min(sm1 - pEff - qPerc - q1, pe*(sm1 < wp ? sm1/wp : 1.0));
		e				= Math.max(0.0, e);
		sm1				= Math.max(0.0, sm1 + p_liquid + sm - pEff - e - qPerc - q1);
		if (sm1 > fc)
		{
			q1			+= sm1 - fc;
			sm1			= fc;
		}
		sm2				= Math.max(0.0, sm2 + qPerc - q2);
		q				= pEff + q1 + q2;
	}
	
	/**
	 * Creates a clone of the element with the current parameter and state value assignments.
	 * Output values are not cloned.
	 * @return A clone of the element with the current parameter and state value assignments
	 */
	public HBVElement cloneState()
	{
		return new HBVElement(tt, dd, fc, beta, pelt, tlt, c, wp, kPerc, k1, k2, s, sm1, sm2);
	}
	
}
