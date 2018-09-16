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
 */

package optimists.vic;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.apache.commons.io.FileUtils;

import maestro_mo.ContVar;
import maestro_mo.MAESTRO;
import maestro_mo.Objective;
import optimists.OPTIMISTS;
import optimists.Particle;
import optimists.ParticleWrapper;
import probDist.KernelDensity;
import probDist.multiVar.KD_GGMLite;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.Point2D;
import utilities.geom.PointSD;
import utilities.stat.ContSeries;
import utilities.thread.Executor;
import utilities.thread.ExecutorThread;
import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumNetwork;
import vic.routing.Simulation;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICAssimilator implements Executor
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public static final String	PARTICLE_ID_PREFIX			= "Particle";
	
	public final static String	OBJ_ID_Q_NSE				= "Streamflow NSE";
	public final static String	OBJ_ID_Q_MAE				= "Streamflow MAE";
	public final static String	OBJ_ID_Q_MARE				= "Streamflow MARE";
	public final static String	OBJ_ID_INDEP_PDF			= "Source mean independent pdf";
	public final static String	OBJ_ID_PDF					= "Source pdf";
	public final static String	OBJ_ID_LOG_PDF				= "Source log-pdf";
	public final static String	OBJ_ID_MAHALANOBIS_DIST		= "Source mean Mahalanobis distance";
	public final static String	OBJ_ID_MAHALANOBIS_FORCE	= "Source Mahalanobis force";
	public final static String	OBJ_ID_MEAN_MAHAL_FORCE		= "Source mean Mahalanobis force";

	public final static String	DATE_TIME_FORMAT_FOLDER		= "yyyy-MM-dd HH.mm";
	public final static String	DATE_TIME_FORMAT_REPORT		= "MM/dd/yyyy HH:mm";
	
	public final static String	REPORT_FOLDER				= "/Reports";
	public final static String	REPORT_FILE 				= "/DA report.txt";
	public final static String	MEAN_Q_FILE					= "/Streamflow.txt";

	private static final String SIMULATION_FAILED 			= "Simulation failed";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------

	private String								parametersFolder;
	private ArrayList<Soil>						soils;
	private ArrayList<String>					globalFileParams;
	private Duration							modelTimeStep;
	private ArrayList<Forcing>					cellForcings;
	private MuskingumNetwork					network;
	private Hashtable<String, Double>			areas;
	private ArrayList<String>					outputs;
	private Hashtable<String, Double>			directFractions;
	
	private String								vicExec;
	private long								simMaxTime;
	private boolean								removeFiles;	
	private Hashtable<LocalDateTime, Double>	obsQ;
	private ArrayList<Objective>				objectives;
	
	private ArrayList<State>					initialStates;
	private String								modelsFolder;
	private NonParametric						currentState;
	
	private HashSet<String>						toDelete;
	private ArrayList<String>					allExecs;
	private ConcurrentLinkedQueue<String>		availExecs;
	
	private String 								meanQFile;
	private LocalDateTime						forecastEnd;
	private String								forecastFolder;
	private ArrayBlockingQueue<ContMultiSample>	forecastQueue;
	private Hashtable<String, String>			forecastQ;
	private Hashtable<LocalDateTime, ContSeries> meanForecastQ;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICAssimilator(String parameterFolder, ArrayList<Soil> soils, 
			ArrayList<String> globalFileParams,	ArrayList<Forcing> cellForcings, 
			Duration modelTimeStep, MuskingumNetwork network, 
			Hashtable<String, Double> areas, ArrayList<String> outputs, 
			Hashtable<String, Double> directFractions, String vicExec, long simMaxTime, 
			boolean removeFiles, Hashtable<LocalDateTime, Double> obsQ, boolean objQNSE, 
			boolean objQMAE, boolean objQMARE, boolean objIndeppdf, boolean objpdf, 
			boolean objLogpdf, boolean objMDist, boolean objMForce, boolean objMeanMForce)
				throws IOException
	{
		this.parametersFolder		= parameterFolder;
		this.soils					= soils;
		this.globalFileParams		= globalFileParams;
		this.cellForcings			= cellForcings;
		this.modelTimeStep			= modelTimeStep;
		this.network				= network;
		this.areas					= areas;
		this.outputs				= outputs;
		this.directFractions		= directFractions;
		
		this.vicExec				= vicExec;
		this.simMaxTime				= simMaxTime;
		this.removeFiles			= removeFiles;
		this.obsQ					= obsQ;
		
		// Create objectives
		objectives					= new ArrayList<>();
		if (objQNSE)
			objectives.add(new Objective(objectives.size(), OBJ_ID_Q_NSE,				true	));
		if (objQMAE)
			objectives.add(new Objective(objectives.size(), OBJ_ID_Q_MAE,				false	));
		if (objQMARE)
			objectives.add(new Objective(objectives.size(), OBJ_ID_Q_MARE,				false	));
		if (objIndeppdf)
			objectives.add(new Objective(objectives.size(), OBJ_ID_INDEP_PDF,			true	));
		if (objpdf)
			objectives.add(new Objective(objectives.size(), OBJ_ID_PDF,					true	));
		if (objLogpdf)
			objectives.add(new Objective(objectives.size(), OBJ_ID_LOG_PDF,				true	));
		if (objMDist)
			objectives.add(new Objective(objectives.size(), OBJ_ID_MAHALANOBIS_DIST,	false	));
		if (objMForce)
			objectives.add(new Objective(objectives.size(), OBJ_ID_MAHALANOBIS_FORCE,	true	));
		if (objMeanMForce)
			objectives.add(new Objective(objectives.size(), OBJ_ID_MEAN_MAHAL_FORCE,	true	));	
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public void assimilate(String problemName, int runIndex, String outputFolder, 
			String modelsFolder, LocalDateTime start, LocalDateTime end, Duration daTimeStep, 
			MAESTRO maestro, ArrayList<State> initialStates, int ensembleSize, int candidateCount,
			int populationSize, double samplePercentage, double rootPercentage, 
			int distType, double scaling, boolean silverman, GGMLiteCreator ggmCreator,
			double particleGreed, int threadCount, boolean assimilatorReports,
			boolean maestroReports, String hallOfFameFolder, long timeLimit) throws IOException
	{
		// Create mean streamflow file
		meanQFile						= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved streamflow (m3/s)\t"
						+ "Mean streamflow (m3/s)\tSt. dev. (m3/s)");
		out.close();
		
		// Create report file folder
		String reportFile				= assimilatorReports ? outputFolder + REPORT_FILE : "";
		String reportFolder				= assimilatorReports ? outputFolder + REPORT_FOLDER : "";
		if (assimilatorReports)
		{
			if (!Files.exists(			FileSystems.getDefault().getPath(reportFolder)))
				Files.createDirectory(	FileSystems.getDefault().getPath(reportFolder));
		}
		
		// Prepare optimizer data
		this.initialStates				= initialStates;
		this.modelsFolder				= modelsFolder;
		ArrayList<ContVar> variables	= initialStates.get(0).getVariables();
		VICParticle defaultParticle		= new VICParticle(this, "default", null, null, null, null,
															null, null);
		String maestroReportFolder		= maestroReports ? outputFolder + "/MAESTRO reports" : "";
		
		// Create initial state distribution
		NonParametric initialState		= null;
		if (distType == OPTIMISTS.TYPE_D_KERNEL || distType == OPTIMISTS.TYPE_F_KERNEL)
		{
			initialState				= new MultiVarKernelDensity();
			initialState.setWeighted(true);
			for (State state : initialStates)
				initialState.addSample(new Sample(1.0, state.toArray()));
			if (distType == OPTIMISTS.TYPE_D_KERNEL)
			{
				if (Double.isNaN(scaling))
					((MultiVarKernelDensity)initialState).computeGaussianDiagBW(silverman);
				else
					((MultiVarKernelDensity)initialState).computeGaussianDiagBW(scaling);
			}	
			else
			{
				if (Double.isNaN(scaling))
					((MultiVarKernelDensity)initialState).computeGaussianBW(silverman);
				else
					((MultiVarKernelDensity)initialState).computeGaussianBW(scaling);
			}
		}
		else if (distType == OPTIMISTS.TYPE_GGM_LITE)
		{
			initialState				= new KD_GGMLite(true);
			for (State state : initialStates)
				initialState.addSample(new Sample(1.0, state.toArray()));
			if (ggmCreator == null)
				((KD_GGMLite)initialState).computeGaussianBW(scaling);
			else
				((KD_GGMLite)initialState).computeGaussianBW(scaling, ggmCreator);
		}
		
		// Create assimilator and parameterize
		OPTIMISTS assimilator			= new OPTIMISTS(problemName, runIndex, defaultParticle, 
												variables, objectives, reportFile, reportFolder, 
												maestroReportFolder, hallOfFameFolder);
		if (maestro != null)
			assimilator.setMaestro(			maestro				);
		assimilator.setEnsembleSize(		ensembleSize		);
		assimilator.setCandidateCount(		candidateCount		);
		assimilator.setPopulationSize(		populationSize		);
		assimilator.setSamplePercentage(	samplePercentage	);
		assimilator.setRootPercentage(		rootPercentage		);
		assimilator.setParticleGreed(		particleGreed		);
		assimilator.setDistributionType(	distType			);
		assimilator.setScaling(				scaling				);
		assimilator.setSilverman(			silverman			);
		assimilator.setThreadCount(			threadCount			);
		
		// Prepare executables
		allExecs				= new ArrayList<>();
		availExecs				= new ConcurrentLinkedQueue<>();
		String execPrefix		= vicExec.substring(0, vicExec.lastIndexOf('.'));
		String extension		= vicExec.substring(vicExec.lastIndexOf('.'), 
													vicExec.length());
		String random			= (int)(10000*Math.random()) + "";
		File source				= new File(vicExec);
		for (int e = 0; e < threadCount; e++)
		{
			String modified		= execPrefix + " " + random + "-" + (e + 1) + extension;
			File target			= new File(modified);
			if (!Files.exists(target.toPath()))
				Files.copy(source.toPath(), target.toPath());
			allExecs.add(modified);
			availExecs.add(modified);
		}
		
		// Run assimilation loop
		LocalDateTime current				= start;
		Duration total						= Duration.between(start, end);
		int totalSteps						= (int) Math.ceil((double)total.toMillis()
														/daTimeStep.toMillis());
		int currentStep						= 0;
		long stepStart						= 0;
		int stepTimeLimit					= Integer.MAX_VALUE;
		long compStart						= System.currentTimeMillis();
		currentState						= initialState;
		DateTimeFormatter formatter			= null;
		toDelete							= new HashSet<>();
		while (!(current.isEqual(end) || current.isAfter(end)))
		{	
			// Create folder
			formatter						= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
			String stepFolder				= modelsFolder + "/" + formatter.format(current);
			if (!Files.exists(FileSystems.getDefault().getPath(stepFolder)))
				Files.createDirectory(FileSystems.getDefault().getPath(stepFolder));
			
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
			currentState		= assimilator.performDATimeStep(current, timeStepEnd, 
																	currentState, stepTimeLimit);
			
			// Store mean streamflow
			Duration daStep_i	= Duration.between(current, timeStepEnd);
			int modelTimeSteps	= (int)(daStep_i.toMinutes()/modelTimeStep.toMinutes());
			writeMeanStreamflow(modelTimeSteps, meanQFile, current, currentState);
			
			// Delete folders
			if (removeFiles)
				try
				{
					FileUtils.deleteDirectory(new File(stepFolder));
				} catch (Exception e)
				{
					toDelete.add(stepFolder);
				}
			
			// Advance control variables
			currentStep++;
			current							= current.plus(daTimeStep);
		}
		
		// Delete remaining folders
		if (removeFiles)
		{
			for (String folder : toDelete)
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e) {}
		}
	}
	
	private void writeMeanStreamflow(int modelTimeSteps, String meanQFile, LocalDateTime current, 
										NonParametric currentState) throws IOException
	{
		PrintWriter out;
		DateTimeFormatter formatter;
		ArrayList<ContSeries> qStats	= new ArrayList<>();
		for (int t = 0; t < modelTimeSteps; t++)
			qStats.add(new ContSeries(true));
		for (int s = 0; s < currentState.getSamples().size(); s++)
		{
			ContMultiSample sample		= currentState.getSamples().get(s);
			ParticleWrapper wrapper		= (ParticleWrapper) sample;
			VICParticle particle		= (VICParticle) wrapper.getParticle();
			ArrayList<Double> q			= particle.getStreamflow();			
			for (int t = 0; t < modelTimeSteps; t++)
			{
				if (q.size() > 0)
				{
					double value		= q.get(t);
					if (!Double.isNaN(value) && !Double.isInfinite(value))
						qStats.get(t).addValue(value, sample.getWeight());
				}
			}
		}
		out 			= new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
		formatter						= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_REPORT);
		LocalDateTime dateTime			= current.plus(modelTimeStep);
		for (int t = 0; t < modelTimeSteps; t++)
		{
			ContSeries stats			= qStats.get(t);
			out.println(formatter.format(dateTime) + "\t" + obsQ.get(dateTime) + "\t" 
										+ stats.getMean() + "\t" + stats.getStDev());
			dateTime					= dateTime.plus(modelTimeStep);
		}
		out.close();
	}
	
	public Particle createNewParticle(int index, LocalDateTime start, LocalDateTime end,
										ArrayList<Double> sourceStateArray)
	{
		// Prepare data
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder				= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		State sourceState				= new State(sourceStateArray, start, initialStates.get(0));
		String error					= "";
		State targetState				= null;
		ArrayList<Double> targetStateArray = null;
		ArrayList<Double> fitness		= null;
		ArrayList<Double> streamF		= null;
		String id						= PARTICLE_ID_PREFIX + " " + index;
		
		// Obtain executable
		String executable				= null;
		while (executable == null)
			executable					= availExecs.poll();
		
		Simulation simulation			= new Simulation(runFolder, sourceState, start, end, 
														modelTimeStep, parametersFolder, soils,
														network, areas,	outputs, directFractions, 
														cellForcings, globalFileParams);
		try
		{
			simulation.run(executable, simMaxTime);
			streamF						= simulation.getStreamflow(outputs.get(0));
			targetState					= simulation.getEndState();
			targetStateArray			= targetState.toArray();
			fitness						= computeFitness(sourceStateArray, streamF, start);
		} catch (Exception e)
		{
			error						= e.getMessage();
			//e.printStackTrace();
			if (error != null)		
				System.out.println("Could not create \"" + id + "\": " + error);
		}
		
		// Stop executable
		String processName				= executable.substring(executable.lastIndexOf('/') + 1, 
															executable.length());
		try
		{
			Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
		} catch (IOException e1)
		{
			e1.printStackTrace();
		}
		availExecs.add(executable);
		
		// Remove simulation files
		if (removeFiles)
		{
			String deleteFolder			= runFolder;
			try
			{
				FileUtils.deleteDirectory(new File(deleteFolder));
				for (String folder : toDelete)
				{
					deleteFolder		= folder;
					FileUtils.deleteDirectory(new File(deleteFolder));
				}
			} catch (Exception e)
			{
				toDelete.add(deleteFolder);
			}
		}
		
		// Verify if there were errors
		if (error == null || !error.equals(""))
		{
			fitness						= new ArrayList<>();
			for (Objective objective : objectives)
				fitness.add(objective.isMaximization()	? Double.NEGATIVE_INFINITY
														: Double.POSITIVE_INFINITY);
		}
		
		// Create particle
		VICParticle particle			= new VICParticle(this, id, sourceStateArray, 
									targetStateArray, sourceState, targetState, streamF, fitness);
		
		// Log
		if (error != null)
		{
			String line					= formatter.format(start) + " - " + id + ": ";
			for (int o = 0; o < objectives.size(); o++)
			{
				Objective objective		= objectives.get(o);
				line					+= objective.getId() + " = ";
				line					+= particle.getFitness(objective.getIndex());
				if (o < objectives.size() - 1)
					line				+= "; ";
			}
			if (forecastFolder == null)
				System.out.println(line);
			if (executable != null)
				availExecs.add(executable);
		}
		
		return particle;
	}	

	private ArrayList<Double> computeFitness(ArrayList<Double> sourceStateArray,
												ArrayList<Double> streamF, LocalDateTime start)
	{
		ArrayList<Double> fitness			= new ArrayList<>();
		for (Objective objective : objectives)
		{
			String objID					= objective.getId();
			boolean ok						= true;
			double fitVal					= Double.NaN;
			double[] source					= Utilities.toArray(sourceStateArray);
			try
			{
				if (objID.equals(OBJ_ID_Q_NSE) || objID.equals(OBJ_ID_Q_MAE) 
						|| objID.equals(OBJ_ID_Q_MARE			))
					fitVal					= getFitnessQ(streamF, start, objID);
				else if (objID.equals(OBJ_ID_INDEP_PDF			))
					fitVal					= currentState.getMeanIndeppdf(				source);
				else if (objID.equals(OBJ_ID_PDF				))
					fitVal					= currentState.getpdf(						source);
				else if (objID.equals(OBJ_ID_LOG_PDF			))
					fitVal					= currentState.getLogpdf(					source);
				else if (objID.equals(OBJ_ID_MAHALANOBIS_DIST	))
					fitVal					= currentState.getMeanMahalanobisDistance(	source);
				else if (objID.equals(OBJ_ID_MAHALANOBIS_FORCE	))
					fitVal					= currentState.getMahalanobisForce(			source);
				else if (objID.equals(OBJ_ID_MEAN_MAHAL_FORCE	))
					fitVal					= currentState.getMeanMahalanobisForce(		source);
				
				if (Double.isNaN(fitVal))
					ok					= false;
				else
					fitness.add(fitVal);
			} catch (Exception e)
			{
				e.printStackTrace();
				ok							= false;
			}
			
			if (!ok)
				fitness.add(objective.isMaximization()	? Double.NEGATIVE_INFINITY
														: Double.POSITIVE_INFINITY);
		}
		return fitness;
	}

	private double getFitnessQ(ArrayList<Double> streamF, LocalDateTime start, String objective)
	{
		// Prepare hydrograph arrays
		ArrayList<Double> modeled		= streamF;
		double[] meaArray				= new double[modeled.size()];
		double[] modArray				= new double[modeled.size()];
		LocalDateTime dateTime			= start.plus(modelTimeStep);
		for (int t = 0; t < modeled.size(); t++)
		{
			meaArray[t]					= obsQ.get(dateTime);
			modArray[t]					= modeled.get(t);
			dateTime					= dateTime.plus(modelTimeStep);
		}
		
		// Compute fitness
		double fitness					= Double.NaN;
		switch (objective)
		{
			case OBJ_ID_Q_NSE:	fitness	= Utilities.computeNashSutcliffe(
												meaArray, modArray);	break;
			case OBJ_ID_Q_MAE:	fitness	= Utilities.computeMeanAbsoluteError(
												meaArray, modArray);	break;
			case OBJ_ID_Q_MARE: fitness	= Utilities.computeMeanAbsRelativeError(
												meaArray, modArray);	break;
		}
		return fitness;
	}

	public String getParticleReportHeader()
	{
		String header	= "";
		for (int o = 0; o < objectives.size(); o++)
		{
			header		+= objectives.get(o).getId();
			if (o < (objectives.size() - 1))
				header	+= "\t";
		}
		return header;
	}

	public String getParticleReport(VICParticle particle)
	{
		String report	= "";
		for (int o = 0; o < objectives.size(); o++)
		{
			report		+= particle.getFitness(o);
			if (o < (objectives.size() - 1))
				report	+= "\t";
		}
		return report;
	}
	
	public void forecast(LocalDateTime forecastEnd, String folder, int threadCount,
			boolean computePerformance) throws IOException
	{
		System.out.println("");
		this.forecastFolder			= folder;
		
		// Prepare mean streamflow
		ParticleWrapper particle	= (ParticleWrapper) currentState.getSamples().get(0);
		LocalDateTime forecastStart	= particle.getEnd().plus(modelTimeStep);
		LocalDateTime current		= forecastStart;
		meanForecastQ				= new Hashtable<>();
		while (!current.isAfter(forecastEnd))
		{
			meanForecastQ.put(current, new ContSeries());
			current					= current.plus(modelTimeStep);
		}
		
		// Populate forecast queue
		int sampleCount				= currentState.getSamples().size();
		Path path					= FileSystems.getDefault().getPath(folder + "/Forecasts");
		if (!Files.exists(path))
			Files.createDirectory(path);
		forecastQueue				= new ArrayBlockingQueue<>(sampleCount);
		for (ContMultiSample sample : currentState.getSamples())
			forecastQueue.offer(sample);
		
		// Launch threads
		this.forecastEnd			= forecastEnd;
		this.forecastQ				= new Hashtable<>();
		for (int t = 0; t < threadCount; t++)
		{
			ExecutorThread thread	= new ExecutorThread(this);
			thread.start("Perform forecast");
		}
		
		// Wait for threads to finish
		long start					= System.currentTimeMillis();
		boolean timeUp				= false;
		while (forecastQ.size() < sampleCount && !timeUp)
			synchronized (this)
			{
				try
				{
					long time		= System.currentTimeMillis() - start;
					timeUp			= time > (sampleCount*simMaxTime)/threadCount;
					wait(50);
				} catch (InterruptedException e) {}
			}
		
		// Write mean streamflow
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_REPORT);
		PrintWriter out 	= new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
		current						= forecastStart;
		Hashtable<LocalDateTime, KernelDensity> dists = new Hashtable<>();
		while (!current.isAfter(forecastEnd))
		{
			ContSeries series		= meanForecastQ.get(current);
			KernelDensity dist = new KernelDensity();
			for (Point2D point : series.getValuesWeights())
				dist.addSample(point.x, point.y);
			dist.computeGaussianBandwidth();
			dists.put(current, dist);
			out.println(formatter.format(current) + "\t" + obsQ.get(current) + "\t" 
									+ series.getMean() + "\t" + series.getStDev());
			current					= current.plus(modelTimeStep);
		}
		out.close();
		
		// Compute forecast performance metrics
		if (computePerformance)
		{
			// Obtain data
			current					= forecastStart;
			ArrayList<Double> obs	= new ArrayList<>();
			ArrayList<Double> mod	= new ArrayList<>();
			ContSeries density		= new ContSeries(false);
			ContSeries rarity		= new ContSeries(false);
			while (!current.isAfter(forecastEnd))
			{
				ContSeries series	= meanForecastQ.get(current);
				Double observed		= obsQ.get(current);
				if (observed != null)
				{
					// Deterministic evaluation
					obs.add(observed);
					mod.add(series.getMean());
					
					// Probabilistic evaluation
					KernelDensity dist = dists.get(current);
					density.addValue(dist.getpdf(observed));
					rarity.addValue(2*Math.abs(dist.getCDF(observed) - 0.5));
				}
				current				= current.plus(modelTimeStep);
			}
			
			// Compute performance
			double[] obsArr			= Utilities.toArray(obs);
			double[] modArr			= Utilities.toArray(mod);
			double nse_l2			= Utilities.computeNashSutcliffe(obsArr, modArr);
			double nse_l1			= Utilities.computeNashSutcliffe(obsArr, modArr, 1.0);
			double mare				= Utilities.computeMeanAbsRelativeError(obsArr, modArr);
			double meanDensity		= density.getMean();
			double meanRarity		= rarity.getMean();
			
			// Store performance
			String file				= folder + "/Forecasts/Performance.txt";
			out					= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
			out.println("Values =\t"			+ obsArr.length	);
			out.println("NSE_l2 =\t"			+ nse_l2		);
			out.println("NSE_l1 =\t"			+ nse_l1		);
			out.println("MARE =\t"				+ mare			);
			out.println("Mean density =\t"		+ meanDensity	);
			out.println("Mean rarity =\t"		+ meanRarity	);
			out.close();
			
			System.out.println("\nValues =\t"		+ obsArr.length	);
			System.out.println("NSE_l2 =\t"			+ nse_l2		);
			System.out.println("NSE_l1 =\t"			+ nse_l1		);
			System.out.println("MARE =\t"			+ mare			);
			System.out.println("Mean density =\t"	+ meanDensity	);
			System.out.println("Mean rarity =\t"	+ meanRarity	);
		}
		
		// Write forecasted streamflow file
		out						= null;
		String file				= folder + "/Forecasts/Streamflow.txt";
		out						= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
		out.println("Model\tWeight\tDischarge (l/s)");
		ArrayList<PointSD> sorter = new ArrayList<>();
		for (ContMultiSample sample : currentState.getSamples())
			sorter.add(new PointSD(((ParticleWrapper) sample).getId(), sample.getWeight(), false));
		Collections.sort(sorter);
		for (int p = sorter.size() - 1; p >= 0; p--)
		{
			String data			= forecastQ.get(sorter.get(p).getX());
			if (data != null)
				out.println(data);
		}
		out.close();
		
		// Remove duplicated executables
		for (String executable : allExecs)
			try
			{
				Files.delete(new File(executable).toPath());
			} catch (IOException e)
			{
				e.printStackTrace();
			}
		
		System.out.println("\nFinished!");
	}
	
	@Override
	public void execute(String processId)
	{
		String executable					= null;
		while (forecastQueue.size() > 0)
		{
			VICParticle particle			= null;
			ContMultiSample sample			= forecastQueue.poll();
			try
			{
				LocalDateTime end			= forecastEnd;
				
				// Prepare model
				ParticleWrapper wrapper		= (ParticleWrapper) sample;
				particle					= (VICParticle) wrapper.getParticle();
				String particleID			= particle.getId();
				LocalDateTime start			= particle.getTargetState().dateTime;
				ArrayList<Double> array		= particle.getTargetStateArray();
				String runFolder			= forecastFolder + "/Forecasts/" + particleID;
				State sourceState			= new State(array, start, initialStates.get(0));
				ArrayList<Double> streamF	= null;
				
				// Obtain executable
				executable					= availExecs.poll();
				
				// Run model
				Simulation simulation		= new Simulation(runFolder, sourceState, start, end, 
										modelTimeStep, parametersFolder, soils, network, areas, 
										outputs, directFractions, cellForcings, globalFileParams);
				simulation.run(executable, simMaxTime);
				
				// Stop executable
				String processName			= executable.substring(executable.lastIndexOf('/') + 1, 
																	executable.length());
				try
				{
					Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
				} catch (IOException e1)
				{
					e1.printStackTrace();
				}
				
				// Retrieve and store hydrograph
				streamF					= simulation.getStreamflow(outputs.get(0));
				double weight			= sample.getWeight();
				String line				= particleID + "\t" + weight + "\t";
				for (Double value : streamF)
					line				+= value + "\t";
				forecastQ.put(particleID, line);
				LocalDateTime current	= start.plus(modelTimeStep);
				int t					= 0;
				while (!current.isAfter(forecastEnd))
				{
					meanForecastQ.get(current).addValue(streamF.get(t), weight);
					current				= current.plus(modelTimeStep);
					t++;
				}
				
				System.out.println("Completed forecast for " + particle.getId());
			} catch (Exception e)
			{
				int attempted			= forecastQ.size();
				String label			= SIMULATION_FAILED + " " + (attempted + 1);
				forecastQ.put(label, label);
				//e.printStackTrace();
				System.out.println("Forecast for " + particle.getId() 
									+ ": " + label + ": " + e.getMessage());
			} finally
			{
				if ( executable != null)
					availExecs.add(executable);
			}
		}
	}
	
}
