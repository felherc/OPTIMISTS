/**
 * This file is part of OPTIMISTS.
 * 
 * OPTIMISTS is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * OPTIMISTS is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with OPTIMISTS. If not,
 * see <https://www.gnu.org/licenses/>.
 * 
 * OPTIMISTS was developed at the University of Pittsburgh by Felipe Hernández under the
 * supervision of Xu Liang, with funding from the U.S. Department of Transportation
 * (award OASRTRS-14-H-PIT) and through the William Kepler Whiteford Professorship of the
 * University of Pittsburgh.
 * 
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 */

package optimists;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.concurrent.ArrayBlockingQueue;

import maestro_mo.MAESTRO;
import maestro_mo.Monitor;
import maestro_mo.Objective;
import maestro_mo.gen.GenWrapper;
import maestro_mo.gen.Generator;
import maestro_mo.pop.Population;
import maestro_mo.pop.groupMerge.Front;
import maestro_mo.pop.groupMerge.GroupMergePopulation;
import maestro_mo.solution.SolutionRoot;
import maestro_mo.solution.SolutionWrapper;
import probDist.Normal;
import probDist.multiVar.EnsembleGGMLite;
import probDist.multiVar.EnsembleNormal;
import probDist.multiVar.KD_GGMLite;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.PointID;

/**
 * Performs data assimilation for a defined time frame from a probability distribution that
 * represents the initial or source state variables and produces a probability distribution for 
 * the final or target state variables
 * 
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class Assimilator implements Monitor
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Default value for {@link #particlesToGenerate}
	 */
	public final static int DEF_PARTICLES_TO_GENERATE = 20;
	
	/**
	 * Default value for {@link #particlesToReturn}
	 */
	public final static int DEF_PARTICLES_TO_RETURN = 20;
	
	/**
	 * Default value for {@link #maxEvaluation}
	 */
	public final static int DEF_MAX_EVALUATION = 1000;
	
	/**
	 * Default value for {@link #samplePercentage}
	 */
	public final static double DEF_SAMPLE_PERCENTAGE = 1.0;
	
	/**
	 * Default value for {@link #rootPercentage}
	 */
	public final static double DEF_ROOT_PERCENTAGE = 0.75;
	
	/**
	 * Default value for {@link #weightingMode}
	 */
	public final static int DEF_WEIGHTING_MODE = OPTIMISTS.WEIGHT_MODE_FRONT;
	
	/**
	 * Default value for {@link #particleGreed}
	 */
	public final static double DEF_PARTICLE_GREED = 0.75;
	
	/**
	 * Label for root particles
	 */
	public static final String LABEL_ROOT = "Root sample";

	/**
	 * Label for randomly sampled particles
	 */
	public static final String LABEL_RANDOM = "Random sample";
	
	/**
	 * Error message: too few samples in source distribution
	 */
	public final static String ERROR_FEW_SAMPLES = "The source distribution should contain at "
														+ "least two samples";
	
	/**
	 * Error message: dimensional mismatch between {@link #sourceStateDist} and {@link #maestro}
	 */
	public final static String ERROR_DIMENSIONAL_MISMATCH = "The sizes of the state distribution "
										+ "and the number of variables in MAESTRO do not match";
	
	/**
	 * Error message: {@link #maestro} is <code>null</code>
	 */
	public final static String ERROR_MAESTRO_NOT_DEFINED = "The MAESTRO optimizer instance was not"
															+ " defined";
	
	/**
	 * Error message: invalid assignment to {@link #weightingMode}
	 */
	public final static String ERROR_INVALID_WEIGHTING_MODE = "Invalid particle weighting mode"
																	+ " selected";

	/**
	 * Wait time while {@link #maestro} is being executed in milliseconds
	 */
	private static final long WAIT_TIME = 1;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Probability distribution of the initial/source state variables
	 */
	private NonParametric sourceStateDist;
	
	/**
	 * The number of particles to generate
	 */
	private int particlesToGenerate;
	
	/**
	 * The number of particles to define the resulting final/target state. Needs to be smaller or
	 * equal to {@link #particlesToGenerate}.
	 */
	private int particlesToReturn;
	
	/**
	 * Maximum number of particles to evaluate (taking into account some might be invalid) and
	 * a different one needs to be created 
	 */
	private int maxEvaluation;
	
	/**
	 * The percentage of {@link #particlesToGenerate} that are to be sampled from the
	 * {@link sourceStateDist}. The complement is generated using {@link #maestro}.
	 */
	private double samplePercentage;
	
	/**
	 * The percentage of the total weight of particles in the {@link sourceStateDist} that is to be
	 * drawn directly. The rest of the particles are sampled randomly from the kernels.
	 */
	private double rootPercentage;
	
	/**
	 * Options available:<ul>
	 * <li>{@link #WEIGHT_MODE_FRONT}
	 * <li>{@link #WEIGHT_MODE_DOMINATION}
	 */
	private int weightingMode;
	
	/**
	 * Controls how much to favor particles in higher-ranked fronts. A value of 0.0 indicates that
	 * all particles are weighted uniformly; a value of 1.0 indicates that only the particles in 
	 * the Pareto Front (those that are non-dominated) are assigned non-zero weights. Conversely, a 
	 * value of -1.0 indicates that only the particles in the last front are assigned non-zero 
	 * weights.
	 */
	private double particleGreed;
	
	/**
	 * The instance of the MAESTRO optimizer to create new particles
	 */
	private MAESTRO maestro;
	
	/**
	 * True if {@link #maestro} finished running and the results can be accessed. False otherwise.
	 */
	private Boolean optimizationCompleted;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param sourceStateDist	{@link #sourceStateDist}
	 */
	public Assimilator(NonParametric sourceStateDist)
	{
		if (sourceStateDist.getSamples().size() < 2)
			throw new IllegalArgumentException(ERROR_FEW_SAMPLES);
		
		this.sourceStateDist		= sourceStateDist;
		this.particlesToGenerate	= DEF_PARTICLES_TO_GENERATE;
		this.particlesToReturn		= DEF_PARTICLES_TO_RETURN;
		this.maxEvaluation			= DEF_MAX_EVALUATION;
		this.samplePercentage		= DEF_SAMPLE_PERCENTAGE;
		this.rootPercentage			= DEF_ROOT_PERCENTAGE;
		this.weightingMode			= DEF_WEIGHTING_MODE;
		this.particleGreed			= DEF_PARTICLE_GREED;
		this.maestro				= null;
		this.optimizationCompleted	= false;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return {@link #particlesToGenerate}
	 */
	public int getParticlesToGenerate()
	{
		return particlesToGenerate;
	}

	/**
	 * @param particlesToGenerate {@link #particlesToGenerate}
	 */
	public void setParticlesToGenerate(int particlesToGenerate)
	{
		this.particlesToGenerate = Math.max(2, particlesToGenerate);
	}

	/**
	 * @return {@link #particlesToReturn}
	 */
	public int getParticlesToReturn()
	{
		return particlesToReturn;
	}

	/**
	 * @param particlesToReturn {@link #particlesToReturn}
	 */
	public void setParticlesToReturn(int particlesToReturn)
	{
		this.particlesToReturn	= Math.max(2, particlesToReturn);
	}

	/**
	 * @return {@link #maxEvaluation}
	 */
	public int getMaxEvaluation()
	{
		return maxEvaluation;
	}

	/**
	 * @param maxEvaluation {@link #maxEvaluation}
	 */
	public void setMaxEvaluation(int maxEvaluation)
	{
		this.maxEvaluation = maxEvaluation;
	}

	/**
	 * @return {@link #samplePercentage}
	 */
	public double getSamplePercentage()
	{
		return samplePercentage;
	}

	/**
	 * @param samplePercentage {@link #samplePercentage}
	 */
	public void setSamplePercentage(double samplePercentage)
	{
		samplePercentage		= Math.max(0.0, samplePercentage);
		this.samplePercentage	= Math.min(1.0, samplePercentage);
	}

	/**
	 * @return {@link #rootPercentage}
	 */
	public double getRootPercentage()
	{
		return rootPercentage;
	}

	/**
	 * @param rootPercentage {@link #rootPercentage}
	 */
	public void setRootPercentage(double rootPercentage)
	{
		rootPercentage		= Math.max(0.0, rootPercentage);
		this.rootPercentage	= Math.min(1.0, rootPercentage);
	}

	/**
	 * @return {@link #weightingMode}
	 */
	public int getWeightingMode()
	{
		return weightingMode;
	}

	/**
	 * @param weightingMode {@link #weightingMode}
	 */
	public void setWeightingMode(int weightingMode)
	{
		this.weightingMode = weightingMode;
	}

	/**
	 * @return {@link #particleGreed}
	 */
	public double getParticleGreed()
	{
		return particleGreed;
	}

	/**
	 * @param particleGreed {@link #particleGreed}
	 */
	public void setParticleGreed(double particleGreed)
	{
		particleGreed		= Math.max(-1.0, particleGreed);
		this.particleGreed	= Math.min( 1.0, particleGreed);
	}
	
	/**
	 * @return {@link #maestro}
	 */
	public MAESTRO getMaestro()
	{
		return maestro;
	}

	/**
	 * @param maestro {@link #maestro}
	 * @throws InvalidParameterException	If the dimensions of {@link #sourceStateDist} and
	 * 										{@link #maestro} do not match
	 */
	public void setMaestro(MAESTRO maestro)
	{
		// Verify values are valid
		if (sourceStateDist.getDimensionality() != maestro.getContVars().size())
			throw new InvalidParameterException(ERROR_DIMENSIONAL_MISMATCH);
		this.maestro = maestro;
	}

	/**
	 * Performs the data assimilation by generating a set of candidate particles, then organizing 
	 * them using non-domination sorting, and then creating a kernel density probability 
	 * distribution that represents the target state variables. The candidate particles are created
	 * by drawing samples from the source probability distribution, randomly generating them from
	 * the distribution, and by using the MAESTRO optimization algorithm.
	 * @param timeLimit		The maximum time to perform the assimilation in milliseconds
	 * @param type			Options for the type of the state probability distribution: 
	 * 						{@link OPTIMISTS#TYPE_D_KERNEL}, {@link OPTIMISTS#TYPE_F_KERNEL}, or
	 * 						{@link OPTIMISTS#TYPE_GGM_LITE}
	 * @param scaling		The factor to scale the sample covariance matrix to produce the 
	 * 						kernels' bandwidth matrix. {@link java.lang.Double#NaN} if either 
	 * 						Silverman's or Scott's rule of thumb should be used instead.
	 * @param silverman		<code>true</code> if the Silverman's rule should be used to determine
	 * 						the scaling factor for the kernels' bandwidth matrix. 
	 * 						<code>false</code> if the Scott's rule should be used instead.
	 * @param dimLimit		The maximum size of the bandwidth matrix (for F-class kernels)
	 * @param corrThreshold	(For F-class kernels) If the absolute value of any computed correlation
	 * 						between state variables is smaller than this threshold, it is assigned
	 * 						a value of zero. That is, use values larger than zero to impose
	 * 						sparsity on the bandwidth matrix.
	 * @param ggmCreator	(For KD-GGMLite) Creator containing the parameters for the GGM kernel
	 * 						of the distribution. <code>null</code> for standard kernels or if 
	 * 						default parameters should be used.
	 * @return The resulting multivariate weighted kernel density probability distribution that
	 * represents the variables of the final/target state
	 */
	public NonParametric assimilate(long timeLimit, int type, double scaling,
				boolean silverman, int dimLimit, double corrThreshold, GGMLiteCreator ggmCreator)
	{
		// Verify optimizer
		if (maestro == null)
			throw new RuntimeException(ERROR_MAESTRO_NOT_DEFINED);
		
		// Draw root samples
		particlesToGenerate					= Math.max(particlesToGenerate, particlesToReturn);
		int sampleCount						= (int) (particlesToGenerate*samplePercentage);
		int totalRoots						= sourceStateDist.getSamples().size();
		ArrayList<ContMultiSample> allRoots	= sourceStateDist.getSamples(true, 1.0, totalRoots);
		int initialRoots					= (int)(sampleCount*rootPercentage);
		initialRoots						= Math.min(allRoots.size(), initialRoots);
		
		// Assign roots to initial sample list and to queue
		ArrayList<ContMultiSample> samples	= new ArrayList<>(particlesToGenerate);
		ArrayBlockingQueue<SolutionRoot> rootQueue
								= new ArrayBlockingQueue<>(Math.max(2, totalRoots - initialRoots));
		for (int r = 0; r < initialRoots; r++)
			samples.add(allRoots.get(r));
		
		// TODO Make this a legitimate option (extra parameter)?
		/*for (int r = samples.size(); r < allRoots.size(); r++)
		{
			ArrayList<Double> values		= allRoots.get(r).getValues();
			rootQueue.add(new SolutionRoot(LABEL_ROOT, null, values, null));
		}*/
		
		// Draw random samples
		int randomCount						= Math.max(0, sampleCount - samples.size());
		double[][] sampleMatrix				= sourceStateDist.sampleMultiple(randomCount);
		for (int s = 0; s < randomCount; s++)
			samples.add(new Sample(1.0, Utilities.toArrayList(sampleMatrix[s])));
		
		// Prepare optimizer
		for (int s = 0; s < samples.size(); s++)
		{
			ContMultiSample sample 	= samples.get(s);
			String label			= s < samples.size() - randomCount ? LABEL_ROOT : LABEL_RANDOM;
			SolutionRoot root		= new SolutionRoot(label, null, sample.getValues(), null);
			maestro.addPredefinedSolution(root);
		}
		ArrayList<GenWrapper> wrappers		= maestro.getGenerators();
		ArrayList<Generator> generators		= new ArrayList<>();
		generators.add(new ParticleGenerator(sourceStateDist, rootQueue));
		for (GenWrapper wrapper : wrappers)
			generators.add(wrapper.getGenerator());
		maestro.setGenerators(generators);
		
		// Run optimization phase
		optimizationCompleted				= false;
		maestro.startOptimization(timeLimit, maxEvaluation, particlesToGenerate);
		while (true)
		{
			// Wait
			try
			{
				synchronized (this)
				{
					wait(WAIT_TIME);
				}
			} catch (Exception e)
			{
				// Do nothing
			}
			
			// Verify if optimization finished
			synchronized (optimizationCompleted)
			{
				if (optimizationCompleted)
					break;
			}
		}
		
		// Verify if the population has enough particles
		Population population				= maestro.getPopulation();
		if (population.size() < particlesToReturn)
		{
			ArrayList<SolutionWrapper> all	= maestro.getAllSolutions();
			if (all != null)
				if (all.size() > population.size())
				{
					ArrayList<Objective> objs = population.getObjectives();
					population				= new GroupMergePopulation();
					population.setCapacity(particlesToReturn);
					for (Objective objective : objs)
						population.addObjective(objective);
					
					// Create copy to avoid thread collisions
					ArrayList<SolutionWrapper> copy = new ArrayList<>();
					if (maestro.getThreadCount() > 1)
						for (int s = 0; s < all.size(); s++)
							copy.add(all.get(s));
					else
						copy				= all;
					
					population.offerSolutions(copy);
					((GroupMergePopulation) population).forceUpdate();
				}
		}
		else
		{
			population.setCapacity(particlesToReturn);
			((GroupMergePopulation) population).forceUpdate();
		}
		
		// Create target distribution
		return computeTargetDistribution((GroupMergePopulation) population, type, scaling,
											silverman, dimLimit, corrThreshold, ggmCreator);
	}
	
	/**
	 * Creates a kernel density probability distribution based on the particles in the provided
	 * population
	 * @param population	The population containing the particles
	 * @param type			Options for the type of the probability distribution: 
	 * 						{@link OPTIMISTS#TYPE_D_KERNEL}, {@link OPTIMISTS#TYPE_F_KERNEL}, or
	 * 						{@link OPTIMISTS#TYPE_GGM_LITE}
	 * @param scaling		The factor to scale the sample covariance matrix to produce the 
	 * 						kernels' bandwidth matrix. {@link java.lang.Double#NaN} if either 
	 * 						Silverman's or Scott's rule of thumb should be used instead.
	 * @param silverman		<code>true</code> if the Silverman's rule should be used to determine
	 * 						the scaling factor for the kernels' bandwidth matrix. 
	 * 						<code>false</code> if the Scott's rule should be used instead.
	 * @param dimLimit		The maximum size of the bandwidth matrix (for F-class kernels)
	 * @param corrThreshold	(For F-class kernels) If the absolute value of any computed correlation
	 * 						is smaller than this threshold, it is assigned a value of zero. That 
	 * 						is, use values larger than zero to impose sparsity on the bandwidth 
	 * 						matrix.
	 * @param ggmCreator	(For KD-GGMLite) Creator containing the parameters for the GGM kernel
	 * 						of the distribution. <code>null</code> for standard kernels or if 
	 * 						default parameters should be used. 
	 * @return The kernel density probability distribution
	 */
	private NonParametric computeTargetDistribution(GroupMergePopulation population, 
			int type, double scaling, boolean silverman, int dimLimit, 
			double corrThreshold, GGMLiteCreator ggmCreator)
	{
		ArrayList<ContMultiSample> samples	= null;
		if (weightingMode == OPTIMISTS.WEIGHT_MODE_FRONT)
			samples							= assignWeightsFront(population, particleGreed);
		else if (weightingMode == OPTIMISTS.WEIGHT_MODE_DOMINATION)
			samples							= assignWeightsDomination(population, particleGreed);
		else
			throw new IllegalArgumentException(ERROR_INVALID_WEIGHTING_MODE);
		
		// Create distribution
		NonParametric dist					= null;
		if (type == OPTIMISTS.TYPE_D_NORMAL || type == OPTIMISTS.TYPE_F_NORMAL)
		{
			dist							= new EnsembleNormal(true);
			dist.setSamples(samples);
			if (type == OPTIMISTS.TYPE_D_NORMAL)
				((EnsembleNormal)dist).computeDiagonalCovariance();
			else
				((EnsembleNormal)dist).computeCovariance(Integer.MAX_VALUE, 0.0);
		}
		else if (type == OPTIMISTS.TYPE_GGM_LITE)
		{
			dist							= new EnsembleGGMLite(true);
			dist.setSamples(samples);
			if (ggmCreator == null)
				((EnsembleGGMLite)dist).computeDependencies();
			else
				((EnsembleGGMLite)dist).computeDependencies(ggmCreator);
		}
		else if (type == OPTIMISTS.TYPE_D_KERNEL || type == OPTIMISTS.TYPE_F_KERNEL)
		{
			dist							= new MultiVarKernelDensity();
			dist.setWeighted(true);
			dist.setSamples(samples);
			if (type == OPTIMISTS.TYPE_D_KERNEL)
			{
				if (Double.isNaN(scaling))
					((MultiVarKernelDensity)dist).computeGaussianDiagBW(silverman);
				else
					((MultiVarKernelDensity)dist).computeGaussianDiagBW(scaling);
			}
			else
			{
				if (Double.isNaN(scaling))
					((MultiVarKernelDensity)dist).computeGaussianBW(silverman, dimLimit,
																	corrThreshold);
				else
					((MultiVarKernelDensity)dist).computeGaussianBW(scaling, dimLimit,
																		corrThreshold);
			}
		}
		else if (type == OPTIMISTS.TYPE_KD_GGM_LITE)
		{
			dist							= new KD_GGMLite(true);
			dist.setSamples(samples);
			if (ggmCreator == null)
				((KD_GGMLite)dist).computeGaussianBW(scaling);
			else
				((KD_GGMLite)dist).computeGaussianBW(scaling, ggmCreator);
		}
		return dist;
	}

	/**
	 * Computes the weights for each of the particles in the provided population according to their
	 * front number
	 * @param population	The population containing the particles in a set of fronts. The
	 * 						weights are computed using the parameters stored within the population
	 * 						to compute from the greed value to the standard deviation of the 
	 * 						sampling Gaussian distribution.
	 * @param greed			A value between -1.0 and 1.0 that represents how much to favor 
	 * 						particles in each of the fronts. 1.0 if only the particles in the first 
	 * 						front are to be assigned non-negligible weights. 0.0 if all particles 
	 * 						should be weighted uniformly. -1.0 if only the particles in the last 
	 * 						front should be assigned non-negligible weights.
	 * @return A list of the weighted samples containing the particles with their computed weight
	 */
	private ArrayList<ContMultiSample> assignWeightsFront(GroupMergePopulation population,
															double greed)
	{
		// Obtain q
		double minQ							= population.getMinQ();
		double maxQ							= population.getMaxQ();
		double greedToQPow					= population.getGreedToQPow();
		greed								= greed >  1.0 ?  1.0 : greed;
		greed								= greed < -1.0 ? -1.0 : greed;
		double temp							= Math.pow(1 - Math.abs(greed), greedToQPow);
		double q							= minQ + (maxQ - minQ)*temp;
		int frontCount						= population.getFronts().size();
		
		// Compute weights
		ArrayList<Double> weights = new ArrayList<>();
		for (int i = 1; i <= frontCount; i++)
			weights.add(Normal.computepdf(1.0, q*frontCount, i));
		
		// Get samples and assign weight
		ArrayList<ContMultiSample> samples	= new ArrayList<>();
		ArrayList<Front> fronts				= population.getFronts();
		for (int f = 0; f < frontCount; f++)
			for (SolutionWrapper sol : fronts.get(f).getSolutions())
			{
				ParticleWrapper particle	= (ParticleWrapper) sol.getSolution();
				double weight				= weights.get(f);
				if (weight > 0.0 && particle.getValues() != null)
				{
					particle.setWeight(weight);
					samples.add(particle);
				}
			}
		return samples;
	}
	
	/**
	 * Computes the weights for each of the particles in the provided population according to the
	 * number of particles each one dominates
	 * @param population	The population containing the particles in a set of fronts
	 * @param greed			A value between -1.0 and 1.0 that represents how much to favor 
	 * 						particles that dominate others. 1.0 if the weight is only defined by
	 * 						the number of particles dominated. 0.0 if all particles should be
	 * 						weighted uniformly. -1.0 the weight should be inversely proportional
	 * 						to the number of particles dominated.
	 * @return A list of the weighted samples containing the particles with their computed weight
	 */
	private ArrayList<ContMultiSample> assignWeightsDomination(GroupMergePopulation population,
														double greed)
	{
		ArrayList<SolutionWrapper> sols	= population.getAllSolutions();
		int solCount					= sols.size();
		ArrayList<ContMultiSample> samples	= new ArrayList<>();
		
		// Uniform weight case
		if (greed == 0.0)
		{
			for (int s = 0; s < solCount; s++)
			{
				ParticleWrapper particle	= (ParticleWrapper) sols.get(s).getSolution();
				if (particle.getValues() != null)
				{
					particle.setWeight(1.0);
					samples.add(particle);
				}
			}
			return samples;
		}
		
		// Prepare domination count
		ArrayList<Integer> dominCount	= new ArrayList<>(solCount);
		for (int s = 0; s < solCount; s++)
			dominCount.add(0);
		ArrayList<Objective> objectives	= population.getObjectives();
		int totalDominationCount		= 0;
		
		// Count times each particle dominates
		for		(int s1 = 0; s1 < solCount - 1; s1++)
		{
			SolutionWrapper sol1		= sols.get(s1);
			for	(int s2 = s1 + 1; s2 < solCount; s2++)
			{
				SolutionWrapper sol2	= sols.get(s2);
				int dominion			= sol1.dominates(sol2, objectives);
				if (dominion > 0)
				{
					int dominationCount	= dominCount.get(s1);
					dominCount.set(s1, dominationCount + 1);
					totalDominationCount++;
				}
				else if (dominion < 0)
				{
					int dominationCount	= dominCount.get(s2);
					dominCount.set(s2, dominationCount + 1);
					totalDominationCount++;
				}
			}
		}
		
		ArrayList<Double> weights		= new ArrayList<>(solCount);
		double baseWeight				= totalDominationCount*(1.0 - Math.abs(greed))/solCount;
		double maxWeight				= Double.NEGATIVE_INFINITY;
		for (int s = 0; s < solCount; s++)
		{
			double weight				= dominCount.get(s);
			weights.add(weight);
			maxWeight					= weight > maxWeight ? weight : maxWeight;
		}
		
		// Get samples and assign weight
		ArrayList<ParticleWrapper> temp	= new ArrayList<>(solCount);
		ArrayList<PointID> order		= new ArrayList<>(solCount);
		for (int s = 0; s < solCount; s++)
		{
			ParticleWrapper particle	= (ParticleWrapper) sols.get(s).getSolution();
			double weight				= weights.get(s);
			if (greed > 0.0)
				weight					= baseWeight + weight;
			if (greed < 0.0)
				weight					= baseWeight + maxWeight - weight;
			if (weight > 0.0 && particle.getValues() != null)
			{
				particle.setWeight(weight);
				temp.add(particle);
				order.add(new PointID(temp.size() - 1, weight, false));
			}
		}
		
		// Sort and return
		Collections.sort(order);
		for (int i = order.size() - 1; i >= 0; i--)
		{
			PointID point				= order.get(i);
			samples.add(temp.get(point.getX()));
		}
		return samples;
	}

	/**
	 * Currently unused.
	 * @return 0
	 */
	public int compare(int objective, ParticleWrapper particleWrapper, ParticleWrapper other)
	{
		return 0;
	}

	@Override
	public void terminate(String message)
	{
		synchronized (optimizationCompleted)
		{
			optimizationCompleted = true;
		}
	}

	@Override
	public void reset()
	{
		// Do nothing
	}

	@Override
	public void populationUpdated()
	{
		// Do nothing
	}
	
}
