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

package optimists.dhsvm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.concurrent.ArrayBlockingQueue;

import org.apache.commons.io.FileUtils;

import dhsvm.MetStation;
import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.Input;
import dhsvm.grid.State;
import dhsvm.stream.StreamNetwork;
import maestro_mo.ContVar;
import maestro_mo.MAESTRO;
import maestro_mo.Objective;
import optimists.OPTIMISTS;
import optimists.ParticleWrapper;
import probDist.KernelDensity;
import probDist.multiVar.EnsembleGGMLite;
import probDist.multiVar.EnsembleNormal;
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

/**
 * Performs multi-objective data assimilation for DHSVM models using OPTIMISTS
 * 
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class DHSVMAssimilator implements Executor
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String	DATE_TIME_FORMAT_FOLDER		= "yyyy-MM-dd HH.mm";
	public final static String	DATE_TIME_FORMAT_FLOW		= "MM.dd.yyyy-HH:mm:ss";
	public final static String	DATE_TIME_FORMAT_REPORT		= "MM/dd/yyyy HH:mm";
	
	public final static String	PARTICLE_ID_PREFIX			= "Particle";
	
	public final static String	OBJ_ID_Q_NSE				= "Streamflow NSE";
	public final static String	OBJ_ID_Q_MAE				= "Streamflow MAE";
	public final static String	OBJ_ID_Q_MARE				= "Streamflow MARE";
	public final static String	OBJ_ID_INDEP_PDF			= "Source mean independent pdf";
	public final static String	OBJ_ID_PDF					= "Source pdf";
	public final static String	OBJ_ID_LOG_PDF				= "Source log-pdf";
	public final static String	OBJ_ID_MAHALANOBIS_DIST		= "Source mean Mahalanobis distance";
	public final static String	OBJ_ID_MAHALANOBIS_FORCE	= "Source Mahalanobis force";
	public final static String	OBJ_ID_MEAN_MAHAL_FORCE		= "Source mean Mahalanobis force";
	public final static String	OBJ_ID_2_TERM_COST			= "2-term cost";
	
	public final static String	REPORT_FOLDER				= "/Reports";
	private static final String REPORT_FILE 				= "/DA report.txt";
	public final static String	MEAN_Q_FILE					= "/Streamflow.txt";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private DHSVMModel							model;
	
	private Duration							modelTimeStep;
	
	private Input								input;
	private int									layers;
	private boolean[][]							osMask;
	private ArrayList<Objective>				objectives;
	private ArrayList<Soil>						soils;
	private StreamNetwork						network;
	
	private ArrayList<String>					metFiles;
	private String								modelsFolder;
	private String								dhsvmExec;
	
	private Hashtable<LocalDateTime, Double>	obsQ;
	private double								obsError;
	private double								bkgrMultiplier;
	
	private NonParametric						currentState;
	
	private String 								meanQFile;
	private LocalDateTime						forecastEnd;
	private String								forecastFolder;
	private ArrayBlockingQueue<ContMultiSample>	forecastQueue;
	private Hashtable<String, String>			forecastQ;
	private Hashtable<LocalDateTime, ContSeries> meanForecastQ;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public DHSVMAssimilator(Duration modelTimeStep, Input input, int layers, ArrayList<Soil> soils,
			ArrayList<Vegetation> vegetations, StreamNetwork network, 
			ArrayList<MetStation> stations, String optionsFile, String areaFile,
			String constantsFile, String dhsvmExec, boolean objQNSE, boolean objQMAE, 
			boolean objQMARE, boolean objIndeppdf, boolean objpdf, boolean objLogpdf,
			boolean objMDist, boolean objMForce, boolean objMeanMForce, boolean obj2TermCost, 
			Hashtable<LocalDateTime, Double> obsQ, double obsError, 
			double bkgrMultiplier) throws IOException
	{
		this.modelTimeStep			= modelTimeStep;
		this.input					= input;
		this.layers					= layers;
		this.soils					= soils;
		this.network				= network;
		this.dhsvmExec				= dhsvmExec;
		this.obsQ					= obsQ;
		this.obsError				= obsError;
		this.bkgrMultiplier			= bkgrMultiplier;
		
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
		if (obj2TermCost)
			objectives.add(new Objective(objectives.size(), OBJ_ID_2_TERM_COST,			false	));
		
		// Copy meteorological files
		metFiles					= new ArrayList<>();
		for (MetStation station : stations)
		{
			String file				= station.dataFile;
			metFiles.add(file);
			station.dataFile		= "../../Met" + file.substring(file.lastIndexOf("/"), 
																		file.length());
		}
		
		// Read configuration sections' files
		String[] optionsSection		= readSectionFile(optionsFile);
		String[] areaSection		= readSectionFile(areaFile);
		String[] constantsSection	= readSectionFile(constantsFile);
		
		// Copy inputs and create model
		String inputFolderCopy		= "../../Input";
		String demFile				= inputFolderCopy + "/DEM.bin";
		String maskFile				= inputFolderCopy + "/Mask.bin";
		String flowDirFile			= inputFolderCopy + "/Flow_dir.bin";
		String soilTypeFile			= inputFolderCopy + "/Soil.bin";
		String soilDepthFile		= inputFolderCopy + "/Soil_depth.bin";
		String vegTypeFile			= inputFolderCopy + "/Vegetation.bin";
		String streamClassFile		= inputFolderCopy + "/stream_class.txt";
		String streamNetworkFile	= inputFolderCopy + "/stream_network.txt";
		String streamMapFile		= inputFolderCopy + "/stream_map.txt";
		String surfaceRoutingFile	= inputFolderCopy + "/surface_routing.txt";
		this.model					= new DHSVMModel(optionsSection, areaSection, 
				constantsSection, demFile, maskFile, flowDirFile, soilTypeFile, soilDepthFile, 
				vegTypeFile, streamClassFile, streamNetworkFile, streamMapFile, 
				surfaceRoutingFile, true, soils, vegetations, stations);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	public void assimilate(String problemName, int runIndex, String outputFolder, 
			String modelsFolder, boolean keepModelFiles, LocalDateTime start, MAESTRO maestro, 
			LocalDateTime end, Duration daTimeStep, ArrayList<State> initialStates, 
			int ensembleSize, int candidateCount, int populationSize, double samplePercentage, 
			double rootPercentage, int distType, double scaling, boolean silverman, 
			GGMLiteCreator ggmCreator, double particleGreed, int threadCount,
			boolean assimilatorReports, boolean maestroReports, String hallOfFameFolder,
			long timeLimit) throws IOException
	{
		// Create model files
		this.modelsFolder				= modelsFolder;
		createModelsFolders(modelsFolder);
		
		// Create mean streamflow file
		meanQFile						= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved streamflow (l/s)\tMean streamflow (l/s)\tSt. dev. (l/s)");
		out.close();
		
		// Create report file folder
		String reportFile				= assimilatorReports ? outputFolder + REPORT_FILE : "";
		String reportFolder				= assimilatorReports ? outputFolder + REPORT_FOLDER : "";
		if (assimilatorReports)
		{
			if (!Files.exists(			FileSystems.getDefault().getPath(reportFolder)))
				Files.createDirectory(	FileSystems.getDefault().getPath(reportFolder));
		}
		
		// Find maximum porosity
		double[] porosity				= new double[3];
		porosity[0]						= Double.NEGATIVE_INFINITY;
		porosity[1]						= Double.NEGATIVE_INFINITY;
		porosity[2]						= Double.NEGATIVE_INFINITY;
		for (Soil soil : soils)
		{
			porosity[0]					= Math.max(soil.getPorosity()[0], porosity[0]);
			porosity[1]					= Math.max(soil.getPorosity()[1], porosity[1]);
			porosity[2]					= Math.max(soil.getPorosity()[2], porosity[2]);
		}
		
		// Obtain over-story mask
		osMask							= new boolean[input.mask.length][input.mask[0].length];
		State sampleState				= initialStates.get(0);
		osMask							= sampleState.getOSMask(input.mask);
		
		// Prepare additional data
		ArrayList<ContVar> variables	= initialStates.get(0).getVariables(input.mask, osMask, 
																				porosity);
		DHSVMParticle defaultParticle	= new DHSVMParticle(this, "default", null, null, null, 
																null, null, null);
		String maestroReportFolder		= maestroReports ? outputFolder + "/MAESTRO reports" : "";
		
		// Create initial state distribution
		NonParametric initialState		= null;
		if (distType == OPTIMISTS.TYPE_D_NORMAL || distType == OPTIMISTS.TYPE_F_NORMAL)
		{
			initialState							= new EnsembleNormal(true);
			for (State state : initialStates)
				initialState.addSample(new Sample(1.0, state.toArray(input.mask, osMask)));
			if (distType == OPTIMISTS.TYPE_D_NORMAL)
				((EnsembleNormal)initialState).computeDiagonalCovariance();
			else
				((EnsembleNormal)initialState).computeCovariance(Integer.MAX_VALUE, 0.0);
		}
		else if (distType == OPTIMISTS.TYPE_GGM_LITE)
		{
			initialState							= new EnsembleGGMLite(true);
			for (State state : initialStates)
				initialState.addSample(new Sample(1.0, state.toArray(input.mask, osMask)));
			if (ggmCreator == null)
				((EnsembleGGMLite)initialState).computeDependencies();
			else
				((EnsembleGGMLite)initialState).computeDependencies(ggmCreator);
		}
		if (distType == OPTIMISTS.TYPE_D_KERNEL || distType == OPTIMISTS.TYPE_F_KERNEL)
		{
			initialState				= new MultiVarKernelDensity();
			initialState.setWeighted(true);
			for (State state : initialStates)
				initialState.addSample(new Sample(1.0, state.toArray(input.mask, osMask)));
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
				initialState.addSample(new Sample(1.0, state.toArray(input.mask, osMask)));
			if (ggmCreator == null)
				((KD_GGMLite)initialState).computeGaussianBW(scaling);
			else
				((KD_GGMLite)initialState).computeGaussianBW(scaling, ggmCreator);
		}
		
		// Create assimilator and parameterize
		OPTIMISTS assimilator	= new OPTIMISTS(problemName, runIndex, defaultParticle, variables, 
					objectives, reportFile, reportFolder, maestroReportFolder, hallOfFameFolder);
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
		ArrayList<String> toDelete			= new ArrayList<>(); 
		while (!(current.isEqual(end) || current.isAfter(end)))
		{			
			// Create folder
			formatter						= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
			String stepFolder				= modelsFolder + "/" + formatter.format(current);
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
			NonParametric newState	= new MultiVarKernelDensity();
			while (newState.getSamples().size() == 0)
				newState					= assimilator.performDATimeStep(current, timeStepEnd, 
												currentState, stepTimeLimit);
			currentState					= newState;
			
			// Store mean streamflow
			writeMeanStreamflow(daTimeStep, meanQFile, current, currentState);
			
			// Delete folders
			if (!keepModelFiles)
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
		
		// Delete folders
		if (!keepModelFiles)
		{
			try
			{
				FileUtils.deleteDirectory(new File(modelsFolder + "/Input"	));
				FileUtils.deleteDirectory(new File(modelsFolder + "/Met"	));
			} catch (Exception e) {}
			
			for (String folder : toDelete)
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e) {}
		}
	}

	private void createModelsFolders(String modelsFolder) throws IOException
	{
		String inputFolder					= modelsFolder + "/Input";
		String metFolder					= modelsFolder + "/Met";
		if (!Files.exists(			FileSystems.getDefault().getPath(inputFolder	)))
			Files.createDirectory(	FileSystems.getDefault().getPath(inputFolder	));
		if (!Files.exists(			FileSystems.getDefault().getPath(metFolder		)))
			Files.createDirectory(	FileSystems.getDefault().getPath(metFolder		));
		input.writeToFiles(inputFolder);
		network.writeToFiles(inputFolder, "stream_class.txt", "stream_network.txt", 
								"stream_map.txt", "surface_routing.txt");
		for (String origFile : metFiles)
		{
			String targetFile				= metFolder + origFile.substring(
												origFile.lastIndexOf("/"), origFile.length());
			File file						= new File(targetFile);
			if (!Files.exists(file.toPath()))
				Files.copy(new File(origFile).toPath(), file.toPath());
		}
	}

	private void writeMeanStreamflow(Duration daTimeStep, String meanQFile, LocalDateTime current, 
										NonParametric currentState) throws IOException
	{
		PrintWriter out;
		DateTimeFormatter formatter;
		ArrayList<ContSeries> qStats	= new ArrayList<>();
		int modelTimeSteps				= (int) daTimeStep.toHours();
		for (int t = 0; t < modelTimeSteps; t++)
			qStats.add(new ContSeries(true));
		for (int s = 0; s < currentState.getSamples().size(); s++)
		{
			ContMultiSample sample		= currentState.getSamples().get(s);
			ParticleWrapper wrapper		= (ParticleWrapper) sample;
			DHSVMParticle particle		= (DHSVMParticle) wrapper.getParticle();
			ArrayList<Double> q			= particle.getStreamflow();
			for (int t = 0; t < modelTimeSteps; t++)
			{
				double value			= q.get(t);
				if (!Double.isNaN(value))
					qStats.get(t).addValue(value, sample.getWeight());
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
	
	public DHSVMParticle createNewParticle(int index, LocalDateTime start, LocalDateTime end,
										ArrayList<Double> sourceStateArray)
	{
		// Prepare data
		State sourceState					= null;
				
		// Prepare files
		DateTimeFormatter formatter			= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder					= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		String stateFolder					= runFolder + "/state";
		String outputFolder					= runFolder + "/output";
		String configFile					= runFolder + "/Configuration.txt";
		try
		{
			if (Files.exists(FileSystems.getDefault().getPath(runFolder	)))
				FileUtils.deleteDirectory(new File(runFolder));
			
			Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
			sourceState						= new State(start, sourceStateArray, input.mask,
												osMask, layers, network.getStreamIDs());
			sourceState.writeToFiles(stateFolder);
			model.writeConfigFile(configFile, start, end, modelTimeStep,
											"state", "output", true, null);
		} catch (IOException e)
		{
			throw new RuntimeException("Could not create model files: " + e.getMessage());
		}
		
		String error						= "";
		State targetState					= null;
		ArrayList<Double> targetStateArray	= null;
		ArrayList<Double> fitness			= null;
		Hashtable<LocalDateTime, Double> modeledQ = null;
		try
		{
			// Run model
			ProcessBuilder pb				= new ProcessBuilder(dhsvmExec, configFile);
			pb.directory(new File(runFolder));
			Process process					= pb.start();
			BufferedReader br				= new BufferedReader(new InputStreamReader(
																	process.getInputStream()));
			while (br.readLine() != null)
			{
				// Do nothing
			}
			
			// Retrieve target state and streamflow
			int rows						= input.elevation.length;
			int cols						= input.elevation[0].length;
			targetState						= new State(end, outputFolder, rows, cols);
			targetStateArray				= targetState.toArray(input.mask, osMask);
			modeledQ 						= getHydrograph(outputFolder);
			
			// Compute fitness values
			fitness							= computeFitness(sourceStateArray, modeledQ);
		} catch (IOException e)
		{
			error							= e.getMessage();
			if (error != null)
				System.out.println("Could not create particle: " + error);
		}
		
		// Verify if there were errors
		if (error == null || !error.equals(""))
		{
			fitness							= new ArrayList<>();
			for (Objective objective : objectives)
				fitness.add(objective.isMaximization()	? Double.NEGATIVE_INFINITY
														: Double.POSITIVE_INFINITY);
		}
		
		// Create particle
		String id							= PARTICLE_ID_PREFIX + " " + index;
		ArrayList<Double> streamFlow		= new ArrayList<>();
		LocalDateTime current				= start.plus(modelTimeStep);
		while (!current.isAfter(end))
		{
			if (modeledQ != null)
				streamFlow.add(modeledQ.get(current));
			else
				streamFlow.add(Double.NaN);
			current							= current.plus(modelTimeStep);
		}
		DHSVMParticle particle				= new DHSVMParticle(this, id, sourceStateArray, 
								targetStateArray, sourceState, targetState, streamFlow, fitness);
		
		// Log
		if (error != null)
		{
			String line						= formatter.format(start) + " - " + id + ": ";
			for (int o = 0; o < objectives.size(); o++)
			{
				Objective objective			= objectives.get(o);
				line						+= objective.getId() + " = ";
				line						+= particle.getFitness(objective.getIndex());
				if (o < objectives.size() - 1)
					line					+= "; ";
			}
			System.out.println(line);
		}
		
		// Return particle
		return particle;
	}

	private ArrayList<Double> computeFitness(ArrayList<Double> sourceStateArray, 
												Hashtable<LocalDateTime, Double> modeledQ)
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
					fitVal					= getFitnessQ(modeledQ, objID);
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
				else if (objID.equals(OBJ_ID_2_TERM_COST))
				{
					double obsCost			= getFitnessQ(modeledQ, objID);
					double[] ssArr			= Utilities.toArray(sourceStateArray);
					double distances[]	 	= currentState.getMahalanobisDistanceToSamples(ssArr);
					double bkgrCost			= 0.0;
					for (int d = 0; d < distances.length; d++)
						bkgrCost			+= distances[d];
					bkgrCost				*= bkgrMultiplier;					
					fitness.add(obsCost + bkgrCost);
				}
				
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
	
	private Hashtable<LocalDateTime, Double> getHydrograph(String outputFolder) 
																	throws FileNotFoundException
	{
		// Obtain modeled hydrograph
		String flowFile					= outputFolder + "/Stream.Flow";
		Scanner scanner					= new Scanner(new FileInputStream(new File(flowFile)));
		Hashtable<LocalDateTime, Double> hydrograph = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double value				= Double.valueOf(tokens[4])/3.6;	// Convert to l/s
			hydrograph.put(dateTime, value);
		}
		scanner.close();
		return hydrograph;
	}

	private double getFitnessQ(Hashtable<LocalDateTime, Double> hydrograph, String objective)
	{		
		// Unpack hydrographs
		ArrayList<Double> modeled		= new ArrayList<>();
		ArrayList<Double> measured		= new ArrayList<>();
		Enumeration<LocalDateTime> keys	= obsQ.size() < hydrograph.size() ? 
											obsQ.keys() : hydrograph.keys();
		while (keys.hasMoreElements())
		{
			LocalDateTime dateTime		= keys.nextElement();
			Double modVal				= hydrograph.get(dateTime);
			Double obsVal				= obsQ.get(dateTime);
			if (modVal != null && obsVal != null)
			{
				modeled.add(modVal);
				measured.add(obsVal);
			}
		}
		
		// Compute fitness
		double[] modArray				= new double[modeled.size()];
		double[] meaArray				= new double[modeled.size()];
		for (int v = 0; v < modeled.size(); v++)
		{
			modArray[v]					= modeled.get(v);
			meaArray[v]					= measured.get(v);
		}
		double fitness					= Double.NaN;
		switch (objective)
		{
			case OBJ_ID_Q_NSE:	fitness	= Utilities.computeNashSutcliffe(
												meaArray, modArray);	break;
			case OBJ_ID_Q_MAE:	fitness	= Utilities.computeMeanAbsoluteError(
												meaArray, modArray);	break;
			case OBJ_ID_Q_MARE: fitness	= Utilities.computeMeanAbsRelativeError(
												meaArray, modArray);	break;
			case OBJ_ID_2_TERM_COST:
				fitness					= 0.0;
				for (int q = 0; q < meaArray.length; q++)
				{
					double diff			= meaArray[q] - modArray[q];
					fitness				+= (1/(obsError*obsError))*diff*diff;
				}
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

	public String getParticleReport(DHSVMParticle particle)
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
	
	private String[] readSectionFile(String textFile) throws FileNotFoundException
	{
		ArrayList<String> lines		= new ArrayList<>();
		Scanner scanner				= new Scanner(new File(textFile));
		while (scanner.hasNextLine())
			lines.add(scanner.nextLine());
		scanner.close();
		
		String[] array				= new String[lines.size()];
		for (int l = 0; l < lines.size(); l++)
			array[l]				= lines.get(l);
		return array;
	}
	
	public void forecast(LocalDateTime forecastEnd, String folder, int threadCount, 
									boolean computePerformance) throws IOException
	{
		System.out.println("");
		this.forecastFolder			= folder;
		if (!folder.equals(modelsFolder))
			createModelsFolders(folder);
		
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
		Files.createDirectory(FileSystems.getDefault().getPath(folder + "/Forecasts"));
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
		while (forecastQ.size() < sampleCount) {}
		
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
			out.println(forecastQ.get(sorter.get(p).getX()));
		out.close();
		out.close();
		System.out.println("\nFinished!");
	}

	@Override
	public void execute(String processID)	// Perform forecast
	{
		while (forecastQueue.size() > 0)
		{
			try
			{
				ContMultiSample sample		= forecastQueue.poll();
				LocalDateTime end			= forecastEnd;
				
				// Prepare model
				ParticleWrapper wrapper		= (ParticleWrapper) sample;
				DHSVMParticle particle		= (DHSVMParticle) wrapper.getParticle();
				String particleID			= particle.getId();
				State finalState			= particle.getTargetState();
				LocalDateTime start			= finalState.dateTime;
				String runFolder			= forecastFolder + "/Forecasts/" + particleID;
				String stateFolder			= runFolder + "/state";
				String outputFolder			= runFolder + "/output";
				String configFile			= runFolder + "/Configuration.txt";
				Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
				
				finalState.writeToFiles(stateFolder);
				model.writeConfigFile(configFile, start, end, modelTimeStep, "state", "output", 
										false, null);
				
				// Run forecast
				ProcessBuilder pb		= new ProcessBuilder(dhsvmExec, configFile);
				pb.directory(new File(runFolder));
				Process process			= pb.start();
				BufferedReader br		= new BufferedReader(new InputStreamReader(
																		process.getInputStream()));
				while (br.readLine() != null)
				{
					// Do nothing
				}
				
				// Retrieve target streamflow
				String line				= particleID + "\t" + sample.getWeight() + "\t";
				Hashtable<LocalDateTime, Double> q = getHydrograph(outputFolder);
				LocalDateTime current 	= start.plus(modelTimeStep);
				while (!current.isAfter(end))
				{
					line				+= q.get(current) + "\t";
					current				= current.plus(modelTimeStep);
				}
				forecastQ.put(particleID, line);
				current					= start.plus(modelTimeStep);
				double weight			= sample.getWeight();
				while (!current.isAfter(forecastEnd))
				{
					meanForecastQ.get(current).addValue(q.get(current), weight);
					current				= current.plus(modelTimeStep);
				}
				System.out.println("Completed forecast for " + particleID);
			} catch (Exception e)
			{
				e.printStackTrace();
			}
		}
	}
	
}
