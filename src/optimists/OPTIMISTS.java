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

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Hashtable;

import maestro_mo.ContVar;
import maestro_mo.MAESTRO;
import maestro_mo.Objective;
import maestro_mo.gen.Generator;
import maestro_mo.solution.SolutionWrapper;
import probDist.multiVar.KD_GGMLite;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;

/**
 * Optimized PareTo Inverse Modeling through Integrated STochastic Search (OPTIMISTS): Performs 
 * data assimilation to compute the final or target probability distribution of a vector 
 * of state variables given an initial or source distribution of the same variables and a 
 * user-provided model. The probability distribution is constructed using a set of candidate 
 * particles that are weighted using multi-objective sorting, and then applying kernel
 * smoothing.
 * 
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class OPTIMISTS
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * State distribution type: diagonal-covariance Normal: {@link EnsembleNormal}
	 */
	public final static int TYPE_D_NORMAL = 0;
	
	/**
	 * State distribution type: full-covariance Normal: {@link EnsembleNormal}
	 */
	public final static int TYPE_F_NORMAL = 1;
	
	/**
	 * State distribution type: Lightweight Gaussian graphical model: {@link EnsembleGGMLite}
	 */
	public final static int TYPE_GGM_LITE = 2;
	
	/**
	 * State distribution type: D-class kernel density: {@link MultiVarKernelDensity} (diagonal 
	 * covariance)
	 */
	public final static int TYPE_D_KERNEL = 3;
	
	/**
	 * State distribution type: F-class kernel density: {@link MultiVarKernelDensity} (full 
	 * covariance)
	 */
	public final static int TYPE_F_KERNEL = 4;
	
	/**
	 * State distribution type: Lightweight Gaussian graphical kernel density: {@link KD_GGMLite}
	 */
	public final static int TYPE_KD_GGM_LITE = 5;
	
	/**
	 * Particle weighting mode in which particles are assigned weights depending on their front in
	 * the population
	 */
	public final static int WEIGHT_MODE_FRONT = 1;
	
	/**
	 * Particle weighting mode in which particles are assigned weights depending on the number of
	 * other particles that each dominates
	 */
	public final static int WEIGHT_MODE_DOMINATION = 2;
	
	/**
	 * Default value for {@link #ensembleSize}
	 */
	public final static int DEF_ENSEMBLE_SIZE = 20;
	
	/**
	 * Default value for {@link #candidateCount}
	 */
	public final static int DEF_CANDIDATE_COUNT = 20;
	
	/**
	 * Default value for {@link #populationSize}
	 */
	public final static int DEF_POPULATION_SIZE = 20;
	
	/**
	 * Default value for {@link #maxEvaluations}
	 */
	public final static int DEF_MAX_EVALUATIONS = 30;
	
	/**
	 * Default value for {@link #samplePercentage}
	 */
	public final static double DEF_SAMPLE_PERCENTAGE = 1.0;
	
	/**
	 * Default value for {@link #rootPercentage}
	 */
	public final static double DEF_ROOT_PERCENTAGE = 0.75;
	
	/**
	 * Default value for {@link #distributionType}
	 */
	public final static int DEF_DISTRIBUTION_TYPE = TYPE_KD_GGM_LITE;
	
	/**
	 * Default value for {@link #scaling}
	 */
	public final static double DEF_SCALING = Double.NaN;
	
	/**
	 * Default value for {@link #silverman}
	 */
	public final static boolean DEF_SILVERMAN = true;
	
	/**
	 * Default value for {@link #dimLimit}
	 */
	public final static int DEF_DIM_LIMIT = Integer.MAX_VALUE;
	
	/**
	 * Default value for {@link #corrThreshold}
	 */
	public final static double DEF_CORR_THRESHOLD = 0.0;
	
	/**
	 * Default value for {@link #particleWeightingMode}
	 */
	public final static int DEF_PARTICLE_WEIGHTING_MODE = WEIGHT_MODE_FRONT;
	
	/**
	 * Default value for {@link #particleGreed}
	 */
	public final static double DEF_PARTICLE_GREED = 0.75;
	
	/**
	 * Default values for {@link #threadCount}
	 */
	public final static int DEF_THREAD_COUNT = 1;
	
	/**
	 * Date time format to name MAESTRO hall of fame files
	 */
	public final static String HALL_OF_FAME_DATE_TIME_FORMAT = "yyyy-MM-dd HH.mm";
	
	/**
	 * Maximum number of variables to be printed in the {@link #reportFile}. If there are more
	 * {@link #stateVariables} than this number, a few statistics are included instead. 
	 */
	public final static int MAX_VARIABLES_TO_PRINT = 200;
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Identifier of the data assimilation problem being addressed
	 */
	private String daProblem;
	
	/**
	 * Index of the current run of the problem
	 */
	private int runIndex;
	
	/**
	 * Implementation of the {@link Particle} interface provided by the user that allows running
	 * the simulations and computing the objective values for candidate states of the model 
	 */
	private Particle defaultParticle;
	
	/**
	 * List with the variables that constitute the states of the model. The order of the variables
	 * is used to organize the parameters of the value list in 
	 * {@link Particle#createNew(int, LocalDateTime, LocalDateTime, ArrayList)}
	 */
	private ArrayList<ContVar> stateVariables;
	
	/**
	 * User-defined list of fitness criteria used to compare candidate state configurations
	 */
	private ArrayList<Objective> objectives;
	
	/**
	 * Instance of the MAESTRO optimization engine containing the customized parameters to use for
	 * the assimilation process, including the selection and parameterization of optimization 
	 * algorithms. <code>null</code> if MAESTRO's default parameters should be used.
	 */
	private MAESTRO maestro;
	
	/**
	 * The number of particles that make up the resulting kernel density probability distributions
	 * of the final/target state variables
	 */
	private int ensembleSize;
	
	/**
	 * The number of candidate particles to generate and evaluate for the creation of 
	 * the resulting probability distributions. Can be the same as {@link #ensembleSize} or larger. 
	 */
	private int candidateCount;
	
	/**
	 * If MAESTRO is to be used to generate candidate particles ({@link #samplePercentage} < 1.0),
	 * indicates the number of candidate particles in the elitist population of the optimizer
	 */
	private int populationSize;
	
	/**
	 * The maximum number of particle evaluations to perform in case some of them fail. Can be the
	 * same as {@link #ensembleSize} or larger.
	 */
	private int maxEvaluations;
	
	/**
	 * The percentage of {@link #candidateCount} that is to be sampled from the initial/source
	 * probability distributions. The complement is generated using MAESTRO.
	 */
	private double samplePercentage;
	
	/**
	 * The percentage of the total weight of particles in the initial/source probability 
	 * distributions that is to be drawn directly. The rest of the particles are sampled randomly 
	 * from the kernels.
	 */
	private double rootPercentage;
	
	/**
	 * Options for the type of the state probability distribution: {@link OPTIMISTS#TYPE_D_KERNEL},
	 * {@link OPTIMISTS#TYPE_F_KERNEL}, or {@link OPTIMISTS#TYPE_KD_GGM_LITE}
	 */
	private int distributionType;
	
	/**
	 * The factor to scale the sample covariance matrix of the particles to produce the bandwidth
	 * matrix of the distribution's kernels. {@link java.lang.Double#NaN} if either Silverman's or
	 * Scott's rule of thumb should be used instead as per {@link #silverman}
	 */
	private double scaling;
	
	/**
	 * <code>true</code> if the Silverman's rule should be used to determine the scaling factor for
	 * the kernels' bandwidth matrix. <code>false</code> if the Scott's rule should be used
	 * instead. Overridden if {@link #scaling} is not {@link java.lang.Double#NaN}.
	 */
	private boolean silverman;
	
	/**
	 * The maximum size of the bandwidth matrix of the kernels (for F-class kernels)
	 */
	private int dimLimit;
	
	/**
	 * If the absolute value of any computed correlation between state variables is smaller than 
	 * this threshold, it is assigned a value of zero. That is, use values larger than zero to 
	 * impose sparsity on the bandwidth matrix of the kernels (for F-class kernels).
	 */
	private double corrThreshold;
	
	/**
	 * Contains the parameters for the GGM kernel of the distribution. <code>null</code> if 
	 * standard kernels are used ({@link #distributionType} = {@link #TYPE_D_KERNEL} or {@link #TYPE_F_KERNEL})
	 * or if default parameters should be selected.
	 */
	private GGMLiteCreator ggmCreator;
	
	/**
	 * Options available:<ul>
	 * <li>{@link #WEIGHT_MODE_FRONT}
	 * <li>{@link #WEIGHT_MODE_DOMINATION}
	 */
	private int particleWeightingMode;
	
	/**
	 * Controls how much to favor particles in higher-ranked fronts. A value of 0.0 indicates that
	 * all particles are weighted uniformly; a value of 1.0 indicates that only the particles in 
	 * the Pareto Front (those that are non-dominated) are assigned non-zero weights. Conversely, a 
	 * value of -1.0 indicates that only the particles in the last front are assigned non-zero 
	 * weights.
	 */
	private double particleGreed;
	
	/**
	 * The number of threads on which to execute the evaluation of the candidate model states
	 * concurrently
	 */
	private int threadCount;
	
	/**
	 * System route to write the general report of the assimilation process. "" if the report 
	 * should not be written.
	 */
	private String reportFile;
	
	/**
	 * System route of the folder to store the reports of each time step of the assimilation 
	 * process. The files include the summary of the particles in the target distribution together
	 * with their fitness values. "" if no such reports should be written.
	 */
	private String reportFolder;
	
	/**
	 * System route of the folder to store the reports of the MAESTRO optimization algorithm of 
	 * each time step of the assimilation process. The files include the configuration of the
	 * optimizer and a list of the solutions in the Pareto front (including their fitness and state 
	 * values). "" if no such reports should be written.
	 */
	private String optimizerReportFolder;
	
	/**
	 * System route of the folder to store the "hall of fame" files for each of the time steps of
	 * the assimilation process. The optimizer updates the file of the current assimilation time
	 * step every time a new solution is promoted to the Pareto front. The file includes all of the
	 * values of the promoted solution. "" if no hall of fame files should be written.
	 */
	private String hallOfFameFolder;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param daProblem				{@link #daProblem}
	 * @param runIndex				{@link #runIndex}
	 * @param defaultParticle		{@link #defaultParticle}
	 * @param stateVariables		{@link #stateVariables}
	 * @param objectives			{@link #objectives}
	 * @param reportFile			{@link #reportFile}
	 * @param reportFolder			{@link #reportFolder}
	 * @param optimizerReportFolder	{@link #optimizerReportFolder}
	 * @param hallOfFameFolder		{@link #hallOfFameFolder}
	 * @throws IOException If the report folders and files cannot be created
	 */
	public OPTIMISTS(String daProblem, int runIndex, Particle defaultParticle, 
						ArrayList<ContVar> stateVariables, ArrayList<Objective> objectives, 
						String reportFile, String reportFolder, String optimizerReportFolder, 
						String hallOfFameFolder) throws IOException
	{
		this.daProblem				= daProblem;
		this.runIndex				= runIndex;
		this.defaultParticle		= defaultParticle;
		this.stateVariables			= stateVariables;
		this.objectives				= objectives;
		this.reportFile				= reportFile;
		this.reportFolder			= reportFolder;
		this.optimizerReportFolder	= optimizerReportFolder;
		this.hallOfFameFolder		= hallOfFameFolder;
		
		// Set defaults
		this.maestro				= null;
		this.ensembleSize			= DEF_ENSEMBLE_SIZE;
		this.candidateCount			= DEF_CANDIDATE_COUNT;
		this.populationSize			= DEF_POPULATION_SIZE;
		this.maxEvaluations			= DEF_MAX_EVALUATIONS;
		this.samplePercentage		= DEF_SAMPLE_PERCENTAGE;
		this.rootPercentage			= DEF_ROOT_PERCENTAGE;
		this.distributionType		= DEF_DISTRIBUTION_TYPE;
		this.scaling				= DEF_SCALING;
		this.silverman				= DEF_SILVERMAN;
		this.dimLimit				= DEF_DIM_LIMIT;
		this.corrThreshold			= DEF_CORR_THRESHOLD;
		this.particleWeightingMode	= DEF_PARTICLE_WEIGHTING_MODE;
		this.particleGreed			= DEF_PARTICLE_GREED;
		this.threadCount			= DEF_THREAD_COUNT;
		
		// Prepare folders
		if (!reportFolder.equals(""))
			if (!Files.exists(FileSystems.getDefault().getPath(			reportFolder			)))
				Files.createDirectory(FileSystems.getDefault().getPath(	reportFolder			));
		if (!optimizerReportFolder.equals(""))
			if (!Files.exists(FileSystems.getDefault().getPath(			optimizerReportFolder	)))
				Files.createDirectory(FileSystems.getDefault().getPath(	optimizerReportFolder	));
		if (!hallOfFameFolder.equals(""))
			if (!Files.exists(FileSystems.getDefault().getPath(			hallOfFameFolder		)))
				Files.createDirectory(FileSystems.getDefault().getPath(	hallOfFameFolder		));
		
		// Write report header
		if (!reportFile.equals(""))
			writeReportHeader();
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods - getters and setters
	// --------------------------------------------------------------------------------------------

	/**
	 * @return {@link #maestro}
	 */
	public MAESTRO getMaestro()
	{
		return maestro;
	}

	/**
	 * @param maestro {@link #maestro}
	 */
	public void setMaestro(MAESTRO maestro)
	{
		this.maestro = maestro;
	}

	/**
	 * @return {@link #ensembleSize}
	 */
	public int getEnsembleSize()
	{
		return ensembleSize;
	}

	/**
	 * @param ensembleSize {@link #ensembleSize}
	 */
	public void setEnsembleSize(int ensembleSize)
	{
		this.ensembleSize = Math.max(2, ensembleSize);
	}

	/**
	 * @return {@link #candidateCount}
	 */
	public int getCandidateCount()
	{
		return candidateCount;
	}

	/**
	 * @param candidateCount {@link #candidateCount}
	 */
	public void setCandidateCount(int candidateCount)
	{
		this.candidateCount = Math.max(2, candidateCount);
	}

	/**
	 * @return {@link #populationSize}
	 */
	public int getPopulationSize()
	{
		return populationSize;
	}

	/**
	 * @param populationSize {@link #populationSize}
	 */
	public void setPopulationSize(int populationSize)
	{
		this.populationSize = Math.max(2, populationSize);
	}
	
	/**
	 * @return {@link #maxEvaluations}
	 */
	public int getMaxEvaluations()
	{
		return maxEvaluations;
	}

	/**
	 * @param maxEvaluations {@link #maxEvaluations}
	 */
	public void setMaxEvaluations(int maxEvaluations)
	{
		this.maxEvaluations = maxEvaluations;
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
	 * @return {@link #distributionType}
	 */
	public int getDistributionType()
	{
		return distributionType;
	}

	/**
	 * @param distributionType {@link #distributionType}
	 */
	public void setDistributionType(int distributionType)
	{
		this.distributionType = distributionType;
	}
	
	/**
	 * @return {@link #scaling}
	 */
	public double getScaling()
	{
		return scaling;
	}

	/**
	 * @param scaling {@link #scaling}
	 */
	public void setScaling(double scaling)
	{
		this.scaling = scaling;
	}

	/**
	 * @return {@link #silverman}
	 */
	public boolean getSilverman()
	{
		return silverman;
	}

	/**
	 * @param silverman {@link #silverman}
	 */
	public void setSilverman(boolean silverman)
	{
		this.silverman = silverman;
	}

	/**
	 * @return {@link #dimLimit}
	 */
	public int getDimLimit()
	{
		return dimLimit;
	}

	/**
	 * @param dimLimit {@link #dimLimit}
	 */
	public void setDimLimit(int dimLimit)
	{
		this.dimLimit = dimLimit;
	}

	/**
	 * @return {@link #corrThreshold}
	 */
	public double getCorrThreshold()
	{
		return corrThreshold;
	}

	/**
	 * @param corrThreshold {@link #corrThreshold}
	 */
	public void setCorrThreshold(double corrThreshold)
	{
		this.corrThreshold = corrThreshold;
	}
	
	/**
	 * @return {@link #particleWeightingMode}
	 */
	public int getParticleWeightingMode()
	{
		return particleWeightingMode;
	}

	/**
	 * @param weightingMode {@link #particleWeightingMode}
	 */
	public void setParticleWeightingMode(int weightingMode)
	{
		this.particleWeightingMode = weightingMode;
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
	 * @return {@link #threadCount}
	 */
	public int getThreadCount()
	{
		return threadCount;
	}

	/**
	 * @param threadCount {@link #threadCount}
	 */
	public void setThreadCount(int threadCount)
	{
		this.threadCount = Math.max(1, threadCount);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	/**
	 * Performs data assimilation for a single time span from the initial/source probability 
	 * distribution of state variables to compute the final/target distribution
	 * @param start			Modeling time stamp of the beginning of the data assimilation time step
	 * @param end			Modeling time stamp of the end of the data assimilation time step
	 * @param sourceState	Probability distribution of the initial or source state variables
	 * @param timeLimit		Time limit for data assimilation in milliseconds
	 * @return				The probability distribution of the final or target state variables
	 */
	public NonParametric performDATimeStep(LocalDateTime start, LocalDateTime end, 
			NonParametric sourceState, int timeLimit)
	{		
		LocalDateTime computationStart		= LocalDateTime.now();
		candidateCount						= Math.max(ensembleSize, candidateCount);
		DateTimeFormatter formatter			= DateTimeFormatter.ISO_LOCAL_DATE_TIME;
		String subProblem					= daProblem + " (" + runIndex + ") - " 
												+ formatter.format(start);
		
		// Create assimilator and particle
		Assimilator assimilator				= new Assimilator(sourceState);
		assimilator.setParticlesToGenerate(	candidateCount			);
		assimilator.setParticlesToReturn(	ensembleSize			);
		assimilator.setMaxEvaluation(		maxEvaluations			);
		assimilator.setSamplePercentage(	samplePercentage		);
		assimilator.setRootPercentage(		rootPercentage			);
		assimilator.setWeightingMode(		particleWeightingMode	);
		assimilator.setParticleGreed(		particleGreed			);
		ParticleWrapper particle			= new ParticleWrapper(defaultParticle, start, end);
		
		// Instantiate MAESTRO
		boolean keepHistory					= populationSize < ensembleSize;
		MAESTRO optimizer					= parameterizeOptimizer(subProblem, particle, 
													assimilator, keepHistory);
		assimilator.setMaestro(optimizer);
		
		// Perform data assimilation
		NonParametric targetState			= assimilator.assimilate(timeLimit,	distributionType,
										scaling, silverman, dimLimit, corrThreshold, ggmCreator);
		
		// Write report files
		try
		{
			if (!reportFile.equals(""))
			{
				PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(
																			reportFile, true)));
				String line		= formatter.format(start) + "\t" + formatter.format(end) + "\t";
				line			+= ensembleSize + "\t" + candidateCount + "\t" + populationSize;
				line			+= "\t" + maxEvaluations + "\t" + samplePercentage + "\t";
				line			+= rootPercentage + "\t" + particleGreed + "\t";
				line			+= formatter.format(computationStart) + "\t" + Duration.between(
											computationStart, LocalDateTime.now()).toMillis();
				out.println(	line);
				out.close();
				
				// TODO Add distribution type and particle weighting mode to report
			}
			if (!reportFolder.equals(""))
				writeStepReport(start, optimizer, targetState);
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		
		// Write MAESTRO report
		try
		{
			if (!optimizerReportFolder.equals(""))
			{
				formatter 		= DateTimeFormatter.ofPattern(HALL_OF_FAME_DATE_TIME_FORMAT);
				String fileName	= optimizerReportFolder + "/" + formatter.format(start) + ".txt";
				System.out.println("Writing optimizer report file " + fileName + "...");
				optimizer.writeReport(true, fileName, true, true, false, false, true);
				System.out.println("Writing file complete!");
			}
		} catch (IOException e)
		{
			e.printStackTrace();
		}
		
		return targetState;
	}
	
	/**
	 * Performs data assimilation from the initial/source probability distribution of state 
	 * variables to compute the final/target distribution. A number of sequential time steps are
	 * used with the provided duration (or less) starting from the start time stamp until the end
	 * time stamp is reached.
	 * @param start			Modeling time stamp of the beginning of the data assimilation time step
	 * @param end			Modeling time stamp of the end of the data assimilation time step
	 * @param daTimeStep	The duration of the assimilation time step (different from the model 
	 * 						time step)
	 * @param sourceState	Probability distribution of the initial or source state variables
	 * @param timeLimit		Time limit for data assimilation in milliseconds
	 * @return				The probability distribution of the final or target state variables
	 */
	public NonParametric performDataAssimilation(LocalDateTime start, LocalDateTime end,
							Duration daTimeStep, NonParametric sourceState, long timeLimit)
	{
		candidateCount						= Math.max(ensembleSize, candidateCount);
		LocalDateTime current				= start;
		Duration total						= Duration.between(start, end);
		int totalSteps						= (int) Math.ceil((double)total.toMillis()
																/daTimeStep.toMillis());
		int currentStep						= 0;
		long stepStart						= 0;
		int stepTimeLimit					= Integer.MAX_VALUE;
		long compStart						= System.currentTimeMillis();
		NonParametric currentState			= sourceState;
		while (!(current.isEqual(end) || current.isAfter(end)))
		{
			// Prepare time step parameters
			LocalDateTime timeStepEnd		= current.plus(daTimeStep);
			timeStepEnd			= timeStepEnd.isAfter(end) ? end : timeStepEnd;
			long now			= System.currentTimeMillis();
			int overTime		= currentStep == 0 ? 0 : (int)(now - stepStart - stepTimeLimit);
			overTime			= Math.max(0, overTime);
			stepStart			= now;
			int timeSoFar		= (int)(stepStart - compStart);
			int timeRemaining	= (int)(timeLimit - timeSoFar);
			stepTimeLimit		= (int)(timeRemaining/(totalSteps - currentStep) - overTime);
			
			// Perform assimilation
			currentState = performDATimeStep(current, timeStepEnd, currentState, stepTimeLimit);
			currentStep++;
			current							= current.plus(daTimeStep);
		}
		return currentState;
	}
	
	/**
	 * Creates an instance of MAESTRO and parameterizes it
	 * @param subProblem	The identifier of the assimilation time step
	 * @param particle		The default particle wrapped with the adequate information
	 * @param assimilator	Instance of the assimilator to assign as the monitor for the optimizer
	 * @param keepHistory	True if the optimizer should maintain all the candidate particles
	 * 						generated
	 * @return The parameterized MAESTRO instance
	 */
	private MAESTRO parameterizeOptimizer(String subProblem, ParticleWrapper particle, 
							Assimilator assimilator, boolean keepHistory)
	{
		MAESTRO optimizer		= new MAESTRO(subProblem, runIndex, particle, assimilator,
												keepHistory, true);
		optimizer.setDiscVars(				null			);
		optimizer.setContVars(				stateVariables	);
		optimizer.setPopulationCapacity(	populationSize	);
		optimizer.setRandomSolutionRatio(	0.0				);
		optimizer.setThreadCount(			threadCount		);
		for (Objective objective : objectives)
		{
			if (objective.isCustom())
				optimizer.addCustomObjective(objective.getIndex(), objective.getId());
			else
				optimizer.addNumericalObjective(objective.getIndex(), objective.getId(),
												objective.isMaximization());						
		}
		if (!hallOfFameFolder.equals(""))
		{
			DateTimeFormatter formatter 	= DateTimeFormatter.ofPattern(
																HALL_OF_FAME_DATE_TIME_FORMAT);
			String hallOfFameFile			= hallOfFameFolder + "/" 
												+ formatter.format(particle.getStart()) + ".txt";
			try
			{
				optimizer.setHallOfFameFile(hallOfFameFile);
			} catch (IOException e)
			{
				e.printStackTrace();
			}
		}
		
		// Assign the parameters from the default instance of MAESTRO
		if (maestro != null)
		{
			optimizer.setConcurrentUpdates(	maestro.getConcurrentUpdates()	);
			optimizer.setGenRatio(			maestro.getGenRatio()			);
			optimizer.setGenMin(			maestro.getGenMin()				);
			optimizer.setAbsGenMin(			maestro.getAbsGenMin()			);
			optimizer.setWeightPop(			maestro.getWeightPop()			);
			optimizer.setWeightFront1(		maestro.getWeightFront1()		);
			for (Generator generator : maestro.getGenerators())
				optimizer.addGenerator(generator);
		}
		
		return optimizer;
	}
	
	/**
	 * Writes the report file for an assimilation time step
	 * @param start			The date and time of the start of the step
	 * @param optimizer		The instance of the optimizer which contains the population with the 
	 * 						best particles
	 * @param targetState	The resulting probability distribution of state variables
	 * @throws IOException	If there is a problem writing the files
	 */
	private void writeStepReport(LocalDateTime start, MAESTRO optimizer, 
									NonParametric targetState) throws IOException
	{
		DateTimeFormatter formatter;
		// Index solutions
		Hashtable<String, SolutionWrapper> index = new Hashtable<>();
		for (int s = 0; s < optimizer.getAllSolutions().size(); s++)
		{
			SolutionWrapper solution = optimizer.getAllSolutions().get(s);
			index.put(solution.getId(), solution);
		}
		
		// Write header
		formatter 			= DateTimeFormatter.ofPattern(HALL_OF_FAME_DATE_TIME_FORMAT);
		String fileName		= reportFolder + "/" + formatter.format(start) + ".txt";
		PrintWriter out 	= new PrintWriter(new BufferedWriter(new FileWriter(fileName, false)));
		String line			= "Particle\tGenerated\tWeight";
		for (Objective objective : objectives)
			line			+= "\t" + objective.getId();
		out.println(line);
		
		// Write particle list
		for (ContMultiSample sample : targetState.getSamples())
		{
			ParticleWrapper wrapper	= (ParticleWrapper) sample;
			SolutionWrapper sol		= index.get(wrapper.getId());
			line					= wrapper.getId() + "\t" + sol.getGeneratorShortId() 
										+ "\t" + sample.getWeight();
			for (Objective objective : objectives)
				line				+= "\t" + wrapper.getFitness(objective.getIndex());
			out.println(line);
		}
		out.close();
	}
	
	/**
	 * Writes the header of the report file
	 * @throws IOException If the file cannot be written
	 */
	private void writeReportHeader() throws IOException
	{
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(reportFile, false)));
		
		// Print file header
		out.println(		"Pareto Data Assimilation execution report"							);
		out.println(		"Software by Felipe Hernández <felher.c@gmail.com>"					);
		out.println(		"Data assimilation problem: " + daProblem							);
		out.println(		"Run index: " + runIndex											);
		out.println(		""																	);
		
		// Print variables
		out.println(		"[State variables]"													);
		if (stateVariables.size() <= MAX_VARIABLES_TO_PRINT)
		{
			out.println(	"Variable\tMinimum\tMaximum"										);
			for (ContVar var : stateVariables)
				out.println(var.getName() + "\t" + var.getMin() + "\t" + var.getMax()			);
		}
		else
			out.println(	"Variable count: " + stateVariables.size()							);
		out.println(		""																	);
		
		// Print objectives
		out.println(		"[Filtering objectives]"											);
		out.println(		"Objective\tCustom?\tTo maximize?"									);
		for (Objective obj : objectives)
			out.println(	obj.getId() + "\t" + obj.isCustom() + "\t" + obj.isMaximization()	);
		out.println(		""																	);
		
		// Print time step table header
		out.println(		"[Assimilation time steps]"											);
		String line			= "Start\tEnd\tEnsemble size\tCandidates\tPopulation size";
		line				+= "\tMax. evaluations\tSample percentage\tRoot percentage";
		line				+= "\\tParticle greed\tComputation start\tTime taken (ms)";
		out.println(		line																);
		
		out.close();
	}
	
}
