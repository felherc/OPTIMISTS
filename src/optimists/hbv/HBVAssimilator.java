package optimists.hbv;

import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Hashtable;

import maestro_mo.ContVar;
import maestro_mo.Objective;
import optimists.OPTIMISTS;
import optimists.Particle;
import probDist.Exponential;
import probDist.KernelDensity;
import probDist.Normal;
import probDist.Uniform;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.Point2D;

public class HBVAssimilator
{

	// --------------------------------------------------------------------------------------------
	// Main
	// --------------------------------------------------------------------------------------------
	
	public static void main(String[] args)
	{
		String problemName			= "Blue River data assimilation";
		int runIndex				= 1;
		
		// Define model parameters
		double tt					= 0.0;
		double dd					= 0.29;
		double fc					= 165.0;
		double beta					= 1.8;
		double pelt					= 7.42;
		double tlt					= -3.14;
		double c					= 0.078;
		double wp					= 22.5;
		double kPerc				= 0.019;
		double k1					= 0.013;
		double k2					= 0.023;		
		double areaMultiplier		= 3031/(24*3.6);
		double s_init				= 50.0;
		double sm1_init				= 2.0;
		double sm2_init				= 7.0;
		HBVElement watershed		= new HBVElement(tt, dd, fc, beta, pelt, tlt, c, wp, kPerc,
														k1, k2, s_init, sm1_init, sm2_init);
		Duration modelTimeStep		= Duration.ofDays(1);
		LocalDateTime daStart		= LocalDateTime.of(1997, 1, 1,  0, 0);
		LocalDateTime daEnd			= LocalDateTime.of(1997, 1, 15, 0, 0);
		
		// Create initial state distribution
		int sampleCount				= 20;
		ArrayList<ContMultiSample> initState = new ArrayList<>(sampleCount);
		Normal snowDist				= new Normal(50.0, 5.0);
		Exponential	sm1Dist			= new Exponential(0.25);
		Uniform		sm2Dist			= new Uniform(6.0, 8.0);
		for (int s = 0; s < sampleCount; s++)
		{
			ArrayList<Double> vals	= new ArrayList<>(3);
			vals.add(snowDist.sample()	);
			vals.add(sm1Dist.sample()	);
			vals.add(sm2Dist.sample()	);
			initState.add(new Sample(1.0, vals));
		}
		
		// Create precipitation and temperature forcings
		Hashtable<LocalDateTime, Point2D> forcing = new Hashtable<>();
		forcing.put(LocalDateTime.of(1997, 1,  1, 0, 0), new Point2D(  0.13,  10.5  ));
		forcing.put(LocalDateTime.of(1997, 1,  2, 0, 0), new Point2D(  0.12,  11.78 ));
		forcing.put(LocalDateTime.of(1997, 1,  3, 0, 0), new Point2D(  0.0,   14.14 ));
		forcing.put(LocalDateTime.of(1997, 1,  4, 0, 0), new Point2D(  0.0,   11.26 ));
		forcing.put(LocalDateTime.of(1997, 1,  5, 0, 0), new Point2D(  0.0,    3.35 ));
		forcing.put(LocalDateTime.of(1997, 1,  6, 0, 0), new Point2D(  0.0,   -0.68 ));
		forcing.put(LocalDateTime.of(1997, 1,  7, 0, 0), new Point2D(  0.06,  -2.21 ));
		forcing.put(LocalDateTime.of(1997, 1,  8, 0, 0), new Point2D(  1.76,  -1.79 ));
		forcing.put(LocalDateTime.of(1997, 1,  9, 0, 0), new Point2D(  0.28,  -2.74 ));
		forcing.put(LocalDateTime.of(1997, 1, 10, 0, 0), new Point2D(  0.0,   -4.65 ));
		forcing.put(LocalDateTime.of(1997, 1, 11, 0, 0), new Point2D(  0.61,  -8.62 ));
		forcing.put(LocalDateTime.of(1997, 1, 12, 0, 0), new Point2D(  0.98, -10.48 ));
		forcing.put(LocalDateTime.of(1997, 1, 13, 0, 0), new Point2D(  0.0,  -10.61 ));
		forcing.put(LocalDateTime.of(1997, 1, 14, 0, 0), new Point2D(  0.0,   -8.78 ));
		
		// Create observed streamflow time series
		Hashtable<LocalDateTime, Double> observedStreamflow = new Hashtable<>();
		observedStreamflow.put(LocalDateTime.of(1997, 1,  2, 0, 0), 7.46);
		observedStreamflow.put(LocalDateTime.of(1997, 1,  3, 0, 0), 7.49); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  4, 0, 0), 7.35); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  5, 0, 0), 7.21); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  6, 0, 0), 6.56); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  7, 0, 0), 6.50); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  8, 0, 0), 6.11); 
		observedStreamflow.put(LocalDateTime.of(1997, 1,  9, 0, 0), 6.45); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 10, 0, 0), 6.62); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 11, 0, 0), 6.67); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 12, 0, 0), 6.59); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 13, 0, 0), 6.22); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 14, 0, 0), 6.28); 
		observedStreamflow.put(LocalDateTime.of(1997, 1, 15, 0, 0), 6.31);
		
		// Define OPTIMISTS parameters
		Duration daTimeStep				= Duration.ofDays(3);
		int threadCount					= 1;
		int ensembleSize				= sampleCount;
		int candidateCount				= sampleCount;
		int populationSize				= 10;
		int maxEvaluations				= 100;
		double samplePercentage			= 1.0;
		double rootPercentage			= 0.6;
		double particleGreed			= 0.5;
		int distType					= OPTIMISTS.TYPE_F_KERNEL;
		double scaling					= 0.1;
		boolean silverman				= true;
		int dimLimit					= Integer.MAX_VALUE;
		boolean weightPerFront			= true;
		
		// Perform assimilation
		HBVAssimilator assimilator = new HBVAssimilator(watershed, areaMultiplier, modelTimeStep);
		ArrayList<ContMultiSample> finalState = assimilator.assimilate(problemName, runIndex,
				daStart, daEnd, daTimeStep, initState, forcing, observedStreamflow,
				ensembleSize, candidateCount, populationSize, maxEvaluations, samplePercentage,
				rootPercentage, dimLimit, distType, scaling, silverman, null, weightPerFront,
				particleGreed, threadCount);
		
		// Create precipitation and temperature forecast forcings
		forcing.put(LocalDateTime.of(1997, 1, 15, 0, 0), new Point2D(  2.08, -5.81 ));
		forcing.put(LocalDateTime.of(1997, 1, 16, 0, 0), new Point2D(  0.0,  -4.36 ));
		forcing.put(LocalDateTime.of(1997, 1, 17, 0, 0), new Point2D(  0.0, -11.54 ));
		forcing.put(LocalDateTime.of(1997, 1, 18, 0, 0), new Point2D(  0.0, -10.46 ));
		forcing.put(LocalDateTime.of(1997, 1, 19, 0, 0), new Point2D(  0.0,  -6.56 ));
		forcing.put(LocalDateTime.of(1997, 1, 20, 0, 0), new Point2D(  0.0,  -0.38 ));
		forcing.put(LocalDateTime.of(1997, 1, 21, 0, 0), new Point2D(  0.43,  6.04 ));
		forcing.put(LocalDateTime.of(1997, 1, 22, 0, 0), new Point2D(  0.82,  6.08 ));
		forcing.put(LocalDateTime.of(1997, 1, 23, 0, 0), new Point2D(  0.0,   4.37 ));
		forcing.put(LocalDateTime.of(1997, 1, 24, 0, 0), new Point2D(  4.52,  5.90 ));
		forcing.put(LocalDateTime.of(1997, 1, 25, 0, 0), new Point2D(  0.0,  -3.11 ));
		forcing.put(LocalDateTime.of(1997, 1, 26, 0, 0), new Point2D(  0.0,  -1.81 ));
		forcing.put(LocalDateTime.of(1997, 1, 27, 0, 0), new Point2D(  0.12, -2.31 ));
		forcing.put(LocalDateTime.of(1997, 1, 28, 0, 0), new Point2D(  0.0, -13.56 ));
		forcing.put(LocalDateTime.of(1997, 1, 29, 0, 0), new Point2D(  0.0, -12.84 ));
		forcing.put(LocalDateTime.of(1997, 1, 30, 0, 0), new Point2D(  0.0,  -7.60 ));
		
		// Perform probabilistic forecast
		LocalDateTime forecastEnd			= LocalDateTime.of(1997, 1, 30, 0, 0);
		ArrayList<KernelDensity> forecast	= assimilator.forecast(daEnd, forecastEnd, finalState,
																	forcing);
		
		// Print forecast results
		System.out.println("\nForecast:");
		System.out.println("Time stamp\tQ (m3/s)\tQ_10%_quantile (m3/s)\tQ_90%_quantile (m3/s)");
		int index							= 0;
		LocalDateTime current				= daEnd;
		while (current.isBefore(forecastEnd))
		{
			KernelDensity dist				= forecast.get(index);
			System.out.println(current.plus(modelTimeStep)
											+ "\t" + dist.getMean()*areaMultiplier
											+ "\t" + dist.getInvCDF(0.1)*areaMultiplier
											+ "\t" + dist.getInvCDF(0.9)*areaMultiplier);
			current							= current.plus(modelTimeStep);
			index++;
		}
	}
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private HBVElement							watershed;
	
	private double								areaMultiplier;
	
	private Duration							modelTimeStep;
	
	private Hashtable<LocalDateTime, Point2D>	forcing;
	
	private Hashtable<LocalDateTime, Double>	observedStreamflow;
	
	private ArrayList<ContVar>					stateVars;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public HBVAssimilator(HBVElement watershed, double areaMultiplier, Duration modelTimeStep)
	{
		this.watershed		= watershed;
		this.areaMultiplier	= areaMultiplier;
		this.modelTimeStep	= modelTimeStep;
		
		// Define state variables
		stateVars	= new ArrayList<>(3);
		stateVars.add(new ContVar("Snow",				0.0,  1000.0			));
		stateVars.add(new ContVar("Soil moisture 1",	0.0, watershed.getFc()	));
		stateVars.add(new ContVar("Soil moisture 2",	0.0, 10000.0			));
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public ArrayList<ContMultiSample> assimilate(String problemName, int runIndex,
			LocalDateTime start, LocalDateTime end, Duration daTimeStep,
			ArrayList<ContMultiSample> initialState, Hashtable<LocalDateTime, Point2D> forcing,
			Hashtable<LocalDateTime, Double> observedStreamflow, int ensembleSize, int candidateCount,
			int populationSize, int maxEvaluations, double samplePercentage, double rootPercentage,
			int dimLimit, int distType, double scaling, boolean silverman,
			GGMLiteCreator ggmCreator, boolean weightPerFront, double particleGreed,
			int threadCount)
	{
		this.forcing					= forcing;
		this.observedStreamflow			= observedStreamflow;

		ArrayList<Objective> objectives	= new ArrayList<>(1);
		objectives.add(new Objective(1, "MAE", false));
		
		HBVParticle defaultParticle		= new HBVParticle(this, "Default", Double.NaN, Double.NaN,
								Double.NaN, Double.NaN, Double.NaN, Double.NaN, null, Double.NaN);
		
		// Create and parameterize assimilator
		OPTIMISTS assimilator			= null;
		try
		{
			assimilator					= new OPTIMISTS(problemName, runIndex, defaultParticle,
															stateVars, objectives, "", "", "", "");
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		int weightingMode				= weightPerFront	? OPTIMISTS.WEIGHT_MODE_FRONT
															: OPTIMISTS.WEIGHT_MODE_DOMINATION;
		assimilator.setEnsembleSize(			ensembleSize		);
		assimilator.setCandidateCount(			candidateCount		);
		assimilator.setPopulationSize(			populationSize		);
		assimilator.setMaxEvaluations(			maxEvaluations		);
		assimilator.setSamplePercentage(		samplePercentage	);
		assimilator.setRootPercentage(			rootPercentage		);
		assimilator.setParticleWeightingMode(	weightingMode		);
		assimilator.setParticleGreed(			particleGreed		);
		assimilator.setDimLimit(				dimLimit			);
		assimilator.setDistributionType(		distType			);
		assimilator.setScaling(					scaling				);
		assimilator.setSilverman(				silverman			);
		assimilator.setThreadCount(				threadCount			);
		
		// Run assimilation loop
		LocalDateTime current			= start;
		ArrayList<ContMultiSample> currentState = initialState;
		while (current.isBefore(end))
		{
			System.out.println("\nPerforming assimilation at " + current.toString() + "...");
			
			LocalDateTime timeStepEnd	= current.plus(daTimeStep);
			currentState				= assimilator.performDATimeStep(current, timeStepEnd,
																currentState, Integer.MAX_VALUE);
			current						= current.plus(daTimeStep);
		}
		
		return currentState;
	}

	public Particle createNewParticle(int index, LocalDateTime start, LocalDateTime end,
			ArrayList<Double> sourceStateArray)
	{
		// Extract values
		String id					= "Particle " + index;
		double snow_init			= sourceStateArray.get(0);
		double soilMoisture1_init	= sourceStateArray.get(1);
		double soilMoisture2_init	= sourceStateArray.get(2);
		
		// Create model
		HBVElement watershedi		= watershed.cloneState();
		watershedi.setS(	snow_init			);
		watershedi.setSm1(	soilMoisture1_init	);
		watershedi.setSm2(	soilMoisture2_init	);
		
		// Run simulation and retrieve hydrographs
		LocalDateTime current		= start;
		ArrayList<Double> modeledQ	= new ArrayList<Double>();
		ArrayList<Double> observedQ	= new ArrayList<Double>();
		while (current.isBefore(end.minus(modelTimeStep)))
		{
			Point2D forcingi		= forcing.get(current);
			double precipitation	= forcingi.getX();
			double temperature		= forcingi.getY();
			
			watershedi.simulateTimeStep(precipitation, temperature);
			modeledQ.add(watershedi.getQ()*areaMultiplier);
			observedQ.add(observedStreamflow.get(current.plus(modelTimeStep)));
			
			current					= current.plus(modelTimeStep);
		}
		
		// Extract final state
		double snow_final			= watershedi.getS();
		double soilMoisture1_final	= watershedi.getSm1();
		double soilMoisture2_final	= watershedi.getSm2();
		
		// Compute objective value
		double[] modArray			= Utilities.toArray(modeledQ);
		double[] obsArray			= Utilities.toArray(observedQ);
		double mae					= Utilities.computeMeanAbsoluteError(obsArray, modArray);
		
		// Print objective and create particle
		System.out.println(id + ": MAE = " + mae);
		return new HBVParticle(this, id, snow_init, soilMoisture1_init, soilMoisture2_init,
						snow_final, soilMoisture1_final, soilMoisture2_final, modArray, mae);
	}
	
	public ArrayList<KernelDensity> forecast(LocalDateTime start, LocalDateTime end,
			ArrayList<ContMultiSample> initialState, Hashtable<LocalDateTime, Point2D> forcing)
	{
		ArrayList<KernelDensity> forecast = new ArrayList<>();
		for (ContMultiSample sample: initialState)
		{
			// Extract initial state
			ArrayList<Double> values	= sample.getValues();
			double snow_init			= values.get(0);
			double soilMoisture1_init	= values.get(1);
			double soilMoisture2_init	= values.get(2);
			
			// Create model
			HBVElement watershedi		= watershed.cloneState();
			watershedi.setS(	snow_init			);
			watershedi.setSm1(	soilMoisture1_init	);
			watershedi.setSm2(	soilMoisture2_init	);
			
			// Run simulation loop
			int index					= 0;
			LocalDateTime current		= start;
			while (current.isBefore(end))
			{
				// Obtain forcing
				Point2D forcingi		= forcing.get(current);
				double precipitation	= forcingi.getX();
				double temperature		= forcingi.getY();
				
				// Run model time step
				watershedi.simulateTimeStep(precipitation, temperature);
				double streamflow		= watershedi.getQ();
				
				// Store streamflow
				if (forecast.size() <= index)
					forecast.add(new KernelDensity());
				KernelDensity qDist		= forecast.get(index);
				qDist.addSample(streamflow, sample.getWeight());
				
				current					= current.plus(modelTimeStep);
				index++;
			}
		}
		
		// Compute bandwidth of distributions and return forecast
		for (KernelDensity dist : forecast)
			dist.computeGaussianBandwidth();
		return forecast;
	}
	
}
