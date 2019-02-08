package optimists.hbv;

import java.time.LocalDateTime;
import java.util.ArrayList;

import optimists.Particle;

/**
 * Represents a candidate particle for an HBV model during a run of OPTIMISTS 
 * @author Felipe Hernández
 */
public class HBVParticle implements Particle
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Instance of the parent assimilator
	 */
	private HBVAssimilator parent;
	
	/**
	 * Identifier
	 */
	private String id;
	
	/**
	 * Initial water-equivalent depth of the snow pack
	 */
	private double snow_init;
	
	/**
	 * Initial water depth in the top soil layer
	 */
	private double soilMoisture1_init;
	
	/**
	 * Initial water depth in the bottom soil layer
	 */
	private double soilMoisture2_init;
	
	/**
	 * Final water-equivalent depth of the snow pack
	 */
	private double snow_final;
	
	/**
	 * Final water depth in the top soil layer
	 */
	private double soilMoisture1_final;
	
	/**
	 * Final water depth in the bottom soil layer
	 */
	private double soilMoisture2_final;
	
	/**
	 * Outflow time series
	 */
	private double[] streamflow;
	
	/**
	 * Mean absolute error after running the model
	 */
	private double mae;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------

	public HBVParticle(HBVAssimilator parent, String id, double snow_init,
			double soilMoisture1_init, double soilMoisture2_init, double snow_final,
			double soilMoisture1_final, double soilMoisture2_final, double[] streamflow,
			double mae)
	{
		this.parent					= parent;
		this.id						= id;
		this.snow_init				= snow_init;
		this.soilMoisture1_init		= soilMoisture1_init;
		this.soilMoisture2_init		= soilMoisture2_init;
		this.snow_final				= snow_final;
		this.soilMoisture1_final	= soilMoisture1_final;
		this.soilMoisture2_final	= soilMoisture2_final;
		this.streamflow				= streamflow;
		this.mae					= mae;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	@Override
	public Particle createNew(int index, LocalDateTime start, LocalDateTime end,
								ArrayList<Double> sourceStateArray)
	{
		return parent.createNewParticle(index, start, end, sourceStateArray);
	}

	@Override
	public String getId()
	{
		return id;
	}

	@Override
	public ArrayList<Double> getSourceStateArray()
	{
		ArrayList<Double> initState = new ArrayList<>(3);
		initState.add(snow_init				);
		initState.add(soilMoisture1_init	);
		initState.add(soilMoisture2_init	);
		return initState;
	}

	@Override
	public ArrayList<Double> getTargetStateArray()
	{
		ArrayList<Double> targetState = new ArrayList<>(3);
		targetState.add(snow_final			);
		targetState.add(soilMoisture1_final	);
		targetState.add(soilMoisture2_final	);
		return targetState;
	}

	@Override
	public String getReportHeader()
	{
		return "MAE";
	}

	@Override
	public String getReport()
	{
		return mae + "";
	}

	@Override
	public double getFitness(int objective)
	{
		return mae;
	}

	@Override
	public int compareTo(int objective, Particle other)
	{
		Double thisMAE	= mae;
		return thisMAE.compareTo(other.getFitness(0));
	}

	public double[] getStreamflow()
	{
		return streamflow;
	}

	@Override
	public boolean isValid()
	{
		return true;
	}
	
}
