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

package optimists.vic.dual;

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
	public final static String	MEAN_VALUES_FILE			= "/Mean values.txt";

	public final static String	SIMULATION_FAILED 			= "Simulation failed";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private String								parametersFolder;
	private ArrayList<String>					globalFileParams;
	private Duration							modelTimeStep;
	private ArrayList<Forcing>					cellForcings;
	private Hashtable<String, Double>			areas;
	private ArrayList<String>					outputs;
	private boolean								defaultParameters;
	
	private String								vicExec;
	private long								simMaxTime;
	private int									maxEvaluationAttempts;
	private int									maxForecastAttempts;
	private boolean								removeFiles;	
	private Hashtable<LocalDateTime, Double>	obsQ;
	private ArrayList<Objective>				objectives;
	
	private ModelConfigurator					configurator;
	
	private String								modelsFolder;
	private NonParametric						currentState;
	private LocalDateTime						simulating;
	
	private HashSet<String>						toDelete;
	private ArrayList<String>					allExecs;
	private ConcurrentLinkedQueue<String>		availExecs;
	
	private String 								meanQFile;
	private LocalDateTime						forecastEnd;
	private String								forecastFolder;
	private ArrayBlockingQueue<ContMultiSample>	forecastQueue;
	private Hashtable<String, String>			forecastQ;
	private Hashtable<String, Integer>			failCount;
	private ArrayList<ContMultiSample>			forecastEndStates;
	
	private Hashtable<LocalDateTime, KernelDensity> distForecastQ;
	private Hashtable<LocalDateTime, KernelDensity> distForecastEv;
	private Hashtable<LocalDateTime, KernelDensity> distForecastSM;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICAssimilator(String parameterFolder, ModelConfigurator configurator,
			ArrayList<String> globalFileParams,	ArrayList<Forcing> cellForcings, 
			Duration modelTimeStep, Hashtable<String, Double> areas, ArrayList<String> outputs,
			boolean defaultParameters, String vicExec, long simMaxTime, int maxEvaluationAttempts,
			int maxForecastAttempts, boolean removeFiles, Hashtable<LocalDateTime,
			Double> obsQ, boolean objQNSE, boolean objQMAE, boolean objQMARE, boolean objIndeppdf,
			boolean objpdf, boolean objLogpdf, boolean objMDist, boolean objMForce,
			boolean objMeanMForce) throws IOException
	{
		this.parametersFolder		= parameterFolder;
		this.configurator			= configurator;
		this.globalFileParams		= globalFileParams;
		this.cellForcings			= cellForcings;
		this.modelTimeStep			= modelTimeStep;
		this.areas					= areas;
		this.outputs				= outputs;
		this.defaultParameters		= defaultParameters;
		
		this.vicExec				= vicExec;
		this.simMaxTime				= simMaxTime;
		this.maxEvaluationAttempts	= maxEvaluationAttempts;
		this.maxForecastAttempts	= maxForecastAttempts;
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

	public NonParametric assimilate(String problemName, int runIndex, String outputFolder, 
			String modelsFolder, LocalDateTime start, LocalDateTime end, Duration daTimeStep, 
			NonParametric initialState, ArrayList<ContVar> variables, MAESTRO maestro,
			int ensembleSize, int candidateCount, int populationSize, int maxEvaluations, 
			double samplePercentage, double rootPercentage, int dimLimit, double corrThreshold,
			int distType, double scaling, boolean silverman, GGMLiteCreator ggmCreator,
			boolean weightPerFront, double particleGreed, int threadCount,
			boolean assimilatorReports, boolean maestroReports, boolean storeTimeSeries,
			String hallOfFameFolder, long timeLimit) throws IOException
	{
		// TODO Remove flag
		System.out.println("\nStarting assimilation...");
		
		// Create mean streamflow file
		meanQFile						= outputFolder + MEAN_VALUES_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved (m3/s)\tMean streamflow (m3/s)\tSt. dev. (m3/s)"
							+ "\t5% quantile\t25% quantile\tmedian\t75% quantile\t95% quantile"
							+ "\tMean evaporation (mm)\tSt. dev. (mm)"
							+ "\t5% quantile\t25% quantile\tmedian\t75% quantile\t95% quantile"
							+ "\tMean soil moisture (mm)\tSt. dev. (mm)"
							+ "\t5% quantile\t25% quantile\tmedian\t75% quantile\t95% quantile");
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
		this.modelsFolder				= modelsFolder;
		VICParticle defaultParticle		= new VICParticle(this, "default", null, null, null, null,
															null, null, null, null, null, false);
		String maestroReportFolder		= maestroReports ? outputFolder + "/MAESTRO reports" : "";
		
		// Create assimilator and parameterize
		int weightingMode				= weightPerFront	? OPTIMISTS.WEIGHT_MODE_FRONT
															: OPTIMISTS.WEIGHT_MODE_DOMINATION;
		OPTIMISTS assimilator			= new OPTIMISTS(problemName, runIndex, defaultParticle, 
												variables, objectives, reportFile, reportFolder, 
												maestroReportFolder, hallOfFameFolder);
		if (maestro != null)
			assimilator.setMaestro(				maestro				);
		assimilator.setEnsembleSize(			ensembleSize		);
		assimilator.setCandidateCount(			candidateCount		);
		assimilator.setPopulationSize(			populationSize		);
		assimilator.setMaxEvaluations(			maxEvaluations		);
		assimilator.setSamplePercentage(		samplePercentage	);
		assimilator.setRootPercentage(			rootPercentage		);
		assimilator.setParticleWeightingMode(	weightingMode		);
		assimilator.setParticleGreed(			particleGreed		);
		assimilator.setThreadCount(				threadCount			);
		assimilator.setDimLimit(				dimLimit			);
		assimilator.setDistributionType(		distType			);
		assimilator.setScaling(					scaling				);
		assimilator.setSilverman(				silverman			);
		assimilator.setCorrThreshold(			corrThreshold		);
		
		// Prepare executables
		if (allExecs == null || allExecs.size() != threadCount)
		{
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
			simulating			= current;
			currentState		= assimilator.performDATimeStep(current, timeStepEnd, 
																	currentState, stepTimeLimit);
			simulating			= null;
			
			// Store mean streamflow
			if (storeTimeSeries)
			{
				Duration daStep_i	= Duration.between(current, timeStepEnd);
				int modelTimeSteps	= (int)(daStep_i.toMinutes()/modelTimeStep.toMinutes());
				writeMeanStreamflow(modelTimeSteps, meanQFile, current, currentState);
			}
			
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
		
		// Store state
		storeState(currentState, variables, outputFolder + "/Final state.txt");
		
		// Delete remaining folders
		if (removeFiles)
		{
			for (String folder : toDelete)
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e) {}
		}
		
		// Remove duplicated executables
		for (String executable : allExecs)
			try
			{
				Files.delete(new File(executable).toPath());
			} catch (IOException e)
			{
				//e.printStackTrace();
			}
		allExecs				= null;
		availExecs				= null;
		
		return currentState;
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
			KernelDensity dist			= new KernelDensity();
			for (Point2D point : stats.getValuesWeights())
				dist.addSample(point.x, point.y);
			dist.computeGaussianBandwidth();
			String line					= formatter.format(dateTime) + "\t" + obsQ.get(dateTime) 
												+ "\t" + stats.getMean() + "\t" + stats.getStDev();
			line						+= "\t" + dist.getInvCDF(0.05) + "\t" 
										+ dist.getInvCDF(0.25) + "\t" + dist.getInvCDF(0.5) + "\t"
										+ dist.getInvCDF(0.75) + "\t" + dist.getInvCDF(0.95);
			out.println(line);
			dateTime					= dateTime.plus(modelTimeStep);
		}
		out.close();
	}
	
	public Particle createNewParticle(int index, LocalDateTime start, LocalDateTime end,
										ArrayList<Double> initialArray)
	{
		// Prepare data
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder				= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		State finalState				= null;
		ArrayList<Double> finalArray	= null;
		ArrayList<Double> fitness		= null;
		ArrayList<Double> streamflow	= null;
		String id						= PARTICLE_ID_PREFIX + " " + index;
		
		// Determine model configuration
		ModelConfiguration conf			= configurator.configure(initialArray, start, 
																	defaultParameters);
		State initialState				= conf.state;
		ArrayList<Soil> soils			= conf.soils;
		MuskingumNetwork network		= conf.network;
		Hashtable<String, Double> directFractions = conf.directFractions;
		Simulation simulation			= new Simulation(runFolder, initialState, start, end,
										modelTimeStep, parametersFolder, soils, network, areas, 
										outputs, directFractions, cellForcings, globalFileParams);
				
		// Prepare simulation attempt loop
		String executable				= availExecs.poll();
		String error					= "";
		boolean done					= false;
		if (executable == null)
		{
			error						= "No more executables available";
			System.out.println(formatter.format(start) + " - " + id + ": " + error);
			done						= true;
		}
		int attempts					= 0;
		while (!done)
		{
			error						= "";
			try
			{
				// Run simulation
				simulation.run(executable, simMaxTime);
				streamflow				= simulation.getStreamflow(outputs.get(0));
				finalState				= simulation.getEndState();
				ModelConfiguration config = new ModelConfiguration(finalState, soils,
																	network, directFractions);
				finalArray				= configurator.toArray(config);
				fitness					= computeFitness(initialArray, streamflow, start);
			} catch (Exception e)
			{
				error					= e.getMessage();
			} finally
			{
				try		// Stop executable
				{
					String processName	= executable.substring(executable.lastIndexOf('/') + 1, 
																executable.length());
					Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
				} catch (Exception e1) {System.out.println("Error: " + e1.getMessage());}
			}
			
			// Verify if the thread should stop
			if (simulating == null)
				error					= "Assimilation finished before completion";
			if (!simulating.isEqual(start))
				error					= "Assimilation restarted before completion";
			
			attempts++;
			if (error.equals(""))
				done					= true;
			else
			{
				done					= attempts >= maxEvaluationAttempts;
				System.out.println(formatter.format(start) + " - " + id + ": " + error
						+ (done ? ", moving on..." : ", re-attempting..."));
				synchronized (this)
				{
					try {wait(200);} catch (InterruptedException e) {}
				}
			}
		}
		
		// Make executable available
		if (executable != null && simulating != null && simulating.isEqual(start))
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
		boolean valid					= true;
		if (error == null || !error.equals(""))
		{
			fitness						= new ArrayList<>();
			for (Objective objective : objectives)
				fitness.add(objective.isMaximization()	? Double.NEGATIVE_INFINITY
														: Double.POSITIVE_INFINITY);
			valid						= false;
		}
		
		// Create particle
		VICParticle particle			= new VICParticle(this, id, initialArray, finalArray,
			initialState, finalState, soils, network, directFractions, streamflow, fitness, valid);
		
		// Log
		if (error.equals(""))
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
			System.out.println(line);
		}
		
		return particle;
	}

	private ArrayList<Double> computeFitness(ArrayList<Double> initialArray, 
												ArrayList<Double> streamflow, LocalDateTime start)
	{
		ArrayList<Double> fitness			= new ArrayList<>();
		for (Objective objective : objectives)
		{
			String objID					= objective.getId();
			boolean ok						= true;
			double fitVal					= Double.NaN;
			double[] source					= Utilities.toArray(initialArray);
			try
			{
				if (objID.equals(OBJ_ID_Q_NSE) || objID.equals(OBJ_ID_Q_MAE) 
						|| objID.equals(OBJ_ID_Q_MARE			))
					fitVal					= getFitnessQ(streamflow, start, objID);
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
	
	public NonParametric getTargetState()
	{
		return currentState;
	}
	
	public void storeState(NonParametric distribution, ArrayList<ContVar> variables, 
							String file) throws IOException
	{
		PrintWriter out			= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
		String header				= "Id\tWeight";
		for (ContVar var : variables)
			header					+= "\t" + var.getName();
		out.println(header);
		for (ContMultiSample sample : distribution.getSamples())
		{
			ParticleWrapper part	= (ParticleWrapper)sample;
			String line				= part.getId() + "\t" + part.getWeight();
			for (Double value : part.getParticle().getTargetStateArray())
				line				+= "\t" + value;
			out.println(line);
		}
		out.close();
	}
	
	public void prepareForecast(NonParametric currentState, LocalDateTime dateTime,
								String outputFolder) throws IOException
	{
		// Create mean streamflow file
		meanQFile							= outputFolder + MEAN_VALUES_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved (m3/s)\tMean streamflow (m3/s)\tSt. dev. (m3/s)"
						+ "\tMean evaporation (mm)\tSt. dev. (mm)"
						+ "\tMean soil moisture (mm)\tSt. dev. (mm)");
		out.close();
		
		// Prepare distribution
		this.currentState					= currentState;
		ArrayList<ContMultiSample> samples	= currentState.getSamples();
		for (int p = 0; p < samples.size(); p++)
		{
			String id						= "Particle " + (p + 1);
			ContMultiSample sample			= samples.get(p);
			ArrayList<Double> valueArray	= sample.getValues();
			ModelConfiguration config		= configurator.configure(valueArray, dateTime,
																		defaultParameters);
			ArrayList<Double> targetArray	= configurator.toArray(config);
			VICParticle particle			= new VICParticle(this, id, null, targetArray, null,
												config.state, config.soils, config.network,
												config.directFractions, null, null, false);
			ParticleWrapper wrapper			= new ParticleWrapper(particle, null, dateTime);
			wrapper.setWeight(sample.getWeight());
			samples.set(p, wrapper);
		}
	}
	
	public NonParametric forecast(LocalDateTime forecastEnd, String folder,
			boolean storeTimeSeries, boolean computePerformance, int threadCount,
			boolean removeFolder) throws IOException
	{
		// TODO Remove flag
		System.out.println("\nStarting forecast...");
		
		this.forecastFolder			= folder;
		
		// Prepare executables
		if (allExecs == null || allExecs.size() != threadCount)
		{
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
		}
		
		// Prepare mean streamflow
		ParticleWrapper particle	= (ParticleWrapper) currentState.getSamples().get(0);
		LocalDateTime forecastStart	= particle.getEnd().plus(modelTimeStep);
		LocalDateTime current		= forecastStart;
		distForecastQ				= new Hashtable<>();
		distForecastEv				= new Hashtable<>();
		distForecastSM				= new Hashtable<>();
		while (!current.isAfter(forecastEnd))
		{
			distForecastQ.put(	current, new KernelDensity());
			distForecastEv.put(	current, new KernelDensity());
			distForecastSM.put(	current, new KernelDensity());
			current					= current.plus(modelTimeStep);
		}
		
		// Populate forecast queue
		int sampleCount				= currentState.getSamples().size();
		forecastEndStates			= new ArrayList<>(sampleCount);
		Path path					= FileSystems.getDefault().getPath(folder + "/Forecasts");
		if (!Files.exists(path))
			Files.createDirectory(path);
		forecastQueue				= new ArrayBlockingQueue<>(sampleCount);
		for (ContMultiSample sample : currentState.getSamples())
			forecastQueue.offer(sample);
		
		// Launch threads
		simulating					= forecastStart.minus(modelTimeStep);
		this.forecastEnd			= forecastEnd;
		this.forecastQ				= new Hashtable<>();
		this.failCount				= new Hashtable<>();
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
					timeUp			= time > (1.5*sampleCount*simMaxTime)/threadCount;
					wait(50);
				} catch (InterruptedException e) {}
			}
		simulating					= null;
		synchronized (this) { try {wait((int)(1.1*simMaxTime)); } catch (InterruptedException e) {} }
		
		// TODO Remove flag
		System.out.println("Finished forecast!");
		
		// Obtain streamflow distributions
		Hashtable<LocalDateTime, KernelDensity> dists = new Hashtable<>();
		current						= forecastStart;
		while (!current.isAfter(forecastEnd))
		{
			KernelDensity distQ		= distForecastQ.get(current);
			distQ.computeGaussianBandwidth();
			dists.put(current, distQ);
			current					= current.plus(modelTimeStep);
		}
		
		// Write mean values
		if (storeTimeSeries)
		{
			DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_REPORT);
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
			current					= forecastStart;
			while (!current.isAfter(forecastEnd))
			{
				// Streamflow
				KernelDensity distQ	= distForecastQ.get(current);
				
				// Evaporation
				KernelDensity distEv = distForecastEv.get(current);
				distEv.computeGaussianBandwidth();
				
				// Soil moisture
				KernelDensity distSM = distForecastSM.get(current);
				distSM.computeGaussianBandwidth();
				
				out.println(formatter.format(current) + "\t" +
								+ distQ.getMean() + "\t" + distQ.getStDev() + "\t" 
								/*+ distQ.getInvCDF(0.05) + "\t" + distQ.getInvCDF(0.25) + "\t"
								+ distQ.getInvCDF(0.5) + "\t" + distQ.getInvCDF(0.75) + "\t"
								+ distQ.getInvCDF(0.95) + "\t"*/
								+ distEv.getMean() + "\t" + distEv.getStDev() + "\t"
								/*+ distEv.getInvCDF(0.05) + "\t" + distEv.getInvCDF(0.25) + "\t"
								+ distEv.getInvCDF(0.5) + "\t" + distEv.getInvCDF(0.75) + "\t"
								+ distEv.getInvCDF(0.95)*/
								+ distSM.getMean() + "\t" + distSM.getStDev() + "\t"
								/*+ distSM.getInvCDF(0.05) + "\t" + distSM.getInvCDF(0.25) + "\t"
								+ distSM.getInvCDF(0.5) + "\t" + distSM.getInvCDF(0.75) + "\t"
								+ distSM.getInvCDF(0.95)*/);
				current				= current.plus(modelTimeStep);
			}
			out.close();
		}
		
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
				KernelDensity dist	= distForecastQ.get(current);
				Double observed		= obsQ.get(current);
				if (observed != null)
				{
					// Deterministic evaluation
					obs.add(observed);
					mod.add(dist.getMean());
					
					// Probabilistic evaluation
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
			PrintWriter out		= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
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
		if (storeTimeSeries)
		{
			String file			= folder + "/Forecasts/Streamflow.txt";
			PrintWriter out		= new PrintWriter(new BufferedWriter(new FileWriter(file, false)));
			out.println("Model\tWeight\tDischarge (l/s)");
			ArrayList<PointSD> sorter = new ArrayList<>();
			for (ContMultiSample sample : currentState.getSamples())
				sorter.add(new PointSD(((ParticleWrapper) sample).getId(),
											sample.getWeight(), false));
			Collections.sort(sorter);
			for (int p = sorter.size() - 1; p >= 0; p--)
			{
				String data		= forecastQ.get(sorter.get(p).getX());
				if (data != null)
					out.println(data);
			}
			out.close();
		}
		
		// Remove duplicated executables
		for (String executable : allExecs)
			try
			{
				Files.delete(new File(executable).toPath());
			} catch (IOException e)
			{
				//e.printStackTrace();
			}
		allExecs				= null;
		availExecs				= null;
	
		// Remove file folder
		if (removeFolder)
		{
			try
			{
				FileUtils.deleteDirectory(new File(path.toString()));
			} catch (Exception e)
			{
				System.out.println("Error removing folder: " + e.getMessage());
			}
		}
		
		// Return final distributions
		return configurator.createDistribution(forecastEndStates);
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getDistForecastQ()
	{
		return distForecastQ;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getDistForecastEv()
	{
		return distForecastEv;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getDistForecastSM()
	{
		return distForecastSM;
	}
	
	public KernelDensity getFinalStreamflow()
	{
		return distForecastQ.get(forecastEnd);
	}
	
	public KernelDensity getFinalEvaporation()
	{
		return distForecastEv.get(forecastEnd);
	}
	
	public KernelDensity getFinalSoilMoisture()
	{
		return distForecastSM.get(forecastEnd);
	}
	
	@Override
	public void execute(String processId)
	{
		while (forecastQueue.size() > 0)
		{
			VICParticle particle			= null;
			String executable				= null;
			LocalDateTime start				= null;
			ContMultiSample sample			= null;
			String particleID				= null; 
			try
			{
				sample						= forecastQueue.poll();
				LocalDateTime end			= forecastEnd;
				
				// Prepare model
				ParticleWrapper wrapper		= (ParticleWrapper) sample;
				particle					= (VICParticle) wrapper.getParticle();
				particleID					= particle.getId();		
				State initialState			= particle.getFinalState();
				ArrayList<Soil> soils		= particle.getSoils();
				MuskingumNetwork network	= particle.getChannels();
				Hashtable<String, Double> directFractions = particle.getDirectFractions();
				start						= initialState.dateTime;
				String runFolder			= forecastFolder + "/Forecasts/" + particleID;
				ArrayList<Double> streamF	= null;
				ArrayList<Double> evap		= null;
				ArrayList<Double> soilM		= null;
				
				// Obtain executable
				executable					= availExecs.poll();
				if (executable == null)
					throw new Exception("No more executables available");
				
				// Run model
				Simulation simulation		= new Simulation(runFolder, initialState, start, end, 
										modelTimeStep, parametersFolder, soils, network, areas, 
										outputs, directFractions, cellForcings, globalFileParams);
				simulation.run(executable, simMaxTime);
				
				// Verify if the thread should stop
				if (simulating == null)
					throw new Exception("Forecast finished before completion");
				if (!simulating.isEqual(start))
					throw new Exception("Forecast restarted before completion");
				
				// Retrieve and store hydrograph
				streamF					= simulation.getStreamflow(outputs.get(0));
				evap					= simulation.getEvaporation();
				soilM					= simulation.getSoilMoisture();
				double weight			= sample.getWeight();				
				String line				= particleID + "\t" + weight + "\t";
				for (Double value : streamF)
					line				+= value + "\t";
				forecastQ.put(particleID, line);
				LocalDateTime current	= start.plus(modelTimeStep);
				int t					= 0;
				while (!current.isAfter(forecastEnd))
				{
					distForecastQ.get(	current).addSample(streamF.get(t),	weight);
					distForecastEv.get(	current).addSample(evap.get(t),		weight);
					distForecastSM.get(	current).addSample(soilM.get(t),	weight);
					current				= current.plus(modelTimeStep);
					t++;
				}
				
				// Retrieve and store final state
				State finalState		= simulation.getEndState();
				ModelConfiguration conf	= new ModelConfiguration(finalState, soils, network,
																	directFractions);
				ArrayList<Double> vals	= configurator.toArray(conf);
				forecastEndStates.add(new Sample(weight, vals));
				
				System.out.println("Completed forecast for " + particleID);
			} catch (Exception e)
			{
				//e.printStackTrace();
				System.out.println("Forecast for " + particleID 
									+ " failed: " + e.getMessage());
				Integer failed			= failCount.get(particleID);
				failed					= 1 + (failed == null ? 0 : failed);
				failCount.put(particleID, failed);
				try
				{
					if (failed >= maxForecastAttempts)
						forecastQ.put(particleID + " attempt " + failed, SIMULATION_FAILED);
					if (simulating != null && simulating.isEqual(start)
							&& failed < maxForecastAttempts)
						forecastQueue.put(sample);
				} catch (InterruptedException e1)
				{
					e1.printStackTrace();
				}
			} finally
			{
				try
				{
					// Stop executable
					String processName	= executable.substring(executable.lastIndexOf('/') + 1, 
																		executable.length());
					Runtime.getRuntime().exec("taskkill /F /IM \"" + processName + "\"");
				} catch (Exception e1)
				{
					//e1.printStackTrace();
					System.out.println("Error stopping executable: " + e1.getMessage());
				}
				
				// Make executable available
				if (executable != null && simulating != null && simulating.isEqual(start))
					availExecs.add(executable);
			}
		}
		
	}
	
}
