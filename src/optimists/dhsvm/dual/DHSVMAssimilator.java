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

package optimists.dhsvm.dual;

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
import java.nio.file.Path;
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
import optimists.Particle;
import optimists.ParticleWrapper;
import optimists.dhsvm.DHSVMModel;
import probDist.KernelDensity;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.PointSD;
import utilities.stat.ContSeries;
import utilities.thread.Executor;
import utilities.thread.ExecutorThread;

/**
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
	public final static String	DATE_TIME_FORMAT_FLOW_1		= "MM.dd.yyyy-HH:mm:ss";
	public final static String	DATE_TIME_FORMAT_FLOW_2		= "MM/dd/yyyy-HH:mm:ss";
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
	
	private DHSVMModel								defaultModel;
	
	private Duration								modelTimeStep;
	
	private Input									input;
	private ArrayList<Objective>					objectives;
	private ArrayList<Soil>							defaultSoils;
	private StreamNetwork							network;
	
	private ModelConfigurator						configurator;
	private boolean									defaultParameters;
	
	private ArrayList<String>						metFiles;
	private String									modelsFolder;
	private String									dhsvmExec;
	
	private Hashtable<LocalDateTime, Double>		obsQ;
	private double									obsError;
	private double									bkgrMultiplier;
	
	private NonParametric							currentState;
	
	private boolean									removeForecastFiles;
	
	private String 									meanQFile;
	private LocalDateTime							forecastEnd;
	private String									forecastFolder;
	private ArrayBlockingQueue<ContMultiSample>		forecastQueue;
	private Hashtable<String, String>				forecastQ;
	private ArrayList<ContMultiSample>				forecastEndStates;
	private ArrayList<LocalDateTime>				statesToSave;
	
	private Hashtable<LocalDateTime, KernelDensity>	forecastStreamflow;
	private Hashtable<LocalDateTime, KernelDensity>	forecastEvaporation;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL1;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL2;
	private Hashtable<LocalDateTime, KernelDensity>	forecastSoilMoistureL3;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public DHSVMAssimilator(ModelConfigurator configurator, boolean defaultParameters, 
			Duration modelTimeStep, Input input, int layers, ArrayList<Soil> soils,
			ArrayList<Vegetation> vegetations, StreamNetwork network,
			ArrayList<MetStation> stations, String optionsFile, String areaFile,
			String constantsFile, String dhsvmExec, boolean objQNSE, boolean objQMAE,
			boolean objQMARE, boolean objIndeppdf, boolean objpdf, boolean objLogpdf,
			boolean objMDist, boolean objMForce, boolean objMeanMForce, boolean obj2TermCost,
			Hashtable<LocalDateTime, Double> obsQ, double obsError, double bkgrMultiplier)
					throws IOException
	{
		this.configurator			= configurator;
		this.defaultParameters		= defaultParameters;
		this.modelTimeStep			= modelTimeStep;
		this.input					= input;
		this.defaultSoils			= soils;
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
		this.defaultModel			= new DHSVMModel(optionsSection, areaSection, 
				constantsSection, demFile, maskFile, flowDirFile, soilTypeFile, soilDepthFile, 
				vegTypeFile, streamClassFile, streamNetworkFile, streamMapFile, 
				surfaceRoutingFile, true, soils, vegetations, stations);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	public NonParametric assimilate(String problemName, int runIndex, String outputFolder, 
			String modelsFolder, boolean keepModelFiles, LocalDateTime start, MAESTRO maestro, 
			LocalDateTime end, Duration daTimeStep, State exampleState, 
			NonParametric initialState, ArrayList<ContVar> variables, int ensembleSize,
			int candidateCount, int populationSize, int maxEvaluations, double samplePercentage,
			double rootPercentage, int dimLimit, int distType, double scaling, boolean silverman,
			GGMLiteCreator ggmCreator, boolean weightPerFront, double particleGreed,
			int threadCount, boolean assimilatorReports, boolean maestroReports,
			String hallOfFameFolder, long timeLimit) throws IOException
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
		for (Soil soil : defaultSoils)
		{
			porosity[0]					= Math.max(soil.getPorosity()[0], porosity[0]);
			porosity[1]					= Math.max(soil.getPorosity()[1], porosity[1]);
			porosity[2]					= Math.max(soil.getPorosity()[2], porosity[2]);
		}
		
		// Prepare additional data
		DHSVMParticle defaultParticle	= new DHSVMParticle(this, "default", null, null, null,
															null, null, null, null, null, null);
		String maestroReportFolder		= maestroReports ? outputFolder + "/MAESTRO reports" : "";
		
		// Create assimilator and parameterize
		int weightingMode				= weightPerFront	? OPTIMISTS.WEIGHT_MODE_FRONT
															: OPTIMISTS.WEIGHT_MODE_DOMINATION;
		OPTIMISTS assimilator	= new OPTIMISTS(problemName, runIndex, defaultParticle, variables, 
					objectives, reportFile, reportFolder, maestroReportFolder, hallOfFameFolder);
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
		assimilator.setDimLimit(				dimLimit			);
		assimilator.setDistributionType(		distType			);
		assimilator.setScaling(					scaling				);
		assimilator.setSilverman(				silverman			);
		assimilator.setThreadCount(				threadCount			);
		
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
		
		System.out.println("");
		
		return currentState;
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

	public Particle createNewParticle(int index, LocalDateTime start, LocalDateTime end,
			ArrayList<Double> valueArray)
	{
		// Prepare data
		State initState						= null;
				
		// Prepare files
		DateTimeFormatter formatter			= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FOLDER);
		String runFolder					= modelsFolder + "/" + formatter.format(start) 
												+ "/Particle " + index;
		String stateFolder					= runFolder + "/state";
		String outputFolder					= runFolder + "/output";
		String configFile					= runFolder + "/Configuration.txt";
		
		// Configure model
		ModelConfiguration config			= configurator.configure(valueArray, start,
																		defaultParameters);
		ArrayList<Double> initArray			= configurator.toArray(config);
		initState							= config.initialState;
		ArrayList<Soil> soils				= config.soils;
		ArrayList<Vegetation> vegetations	= config.vegetations;
		StreamNetwork network				= config.network;
		String classFileName				= "stream_class.txt";
		String networkFileName				= "stream_network.txt";
		String mapFileName					= "stream_map.txt";
		String surfaceRoutingFile			= "surface_routing.txt";
		DHSVMModel model					= new DHSVMModel(defaultModel.optionsSection,
				defaultModel.areaSection, defaultModel.constantsSection, defaultModel.demFile,
				defaultModel.maskFile, defaultModel.flowDirFile, defaultModel.soilTypeFile,
				defaultModel.soilDepthFile, defaultModel.vegTypeFile, classFileName,
				networkFileName, mapFileName, surfaceRoutingFile, defaultModel.d8FlowDirection,
				soils, vegetations, defaultModel.stations);
		
		// Write files
		try
		{
			if (Files.exists(FileSystems.getDefault().getPath(runFolder	)))
				FileUtils.deleteDirectory(new File(runFolder));
			
			Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
			Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
			network.writeToFiles(runFolder, classFileName, networkFileName, mapFileName,
									surfaceRoutingFile);
			
			initState.writeToFiles(stateFolder);
			model.writeConfigFile(configFile, start, end, modelTimeStep, "state", "output", true, null);
		} catch (IOException e)
		{
			throw new RuntimeException("Could not create model files: " + e.getMessage());
		}
		
		String error						= "";
		State endState						= null;
		ArrayList<Double> endArray			= null;
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
				synchronized (this)	{try {wait(1);} catch (InterruptedException e) {}}
			}
			
			// Retrieve target state and streamflow
			int rows						= input.elevation.length;
			int cols						= input.elevation[0].length;
			endState						= new State(end, outputFolder, rows, cols);
			endArray						= configurator.toArray(new ModelConfiguration(
														soils, vegetations, network, endState));
			modeledQ 						= getHydrograph(outputFolder);
			
			// Compute fitness values
			fitness							= computeFitness(valueArray, modeledQ);
		} catch (IOException e)
		{
			String id						= PARTICLE_ID_PREFIX + " " + index;
			String line						= formatter.format(start) + " - " + id + ": ";
			error							= e.getMessage();
			line							+= error;
			System.out.println(line);
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
		ArrayList<Double> streamflow		= new ArrayList<>();
		LocalDateTime current				= start.plus(modelTimeStep);
		while (!current.isAfter(end))
		{
			if (modeledQ != null)
				streamflow.add(modeledQ.get(current));
			else
				streamflow.add(Double.NaN);
			current							= current.plus(modelTimeStep);
		}
		DHSVMParticle particle				= new DHSVMParticle(this, id, soils, vegetations,
						network, initArray, endArray, initState, endState, streamflow, fitness);
		
		// Log
		if (error != null && error.equals(""))
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
	
	private Hashtable<LocalDateTime, Double> getHydrograph(String outputFolder) 
																	throws FileNotFoundException
	{
		// Obtain modeled hydrograph
		String flowFile					= outputFolder + "/Stream.Flow";
		Scanner scanner					= new Scanner(new FileInputStream(new File(flowFile)));
		Hashtable<LocalDateTime, Double> hydrograph = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_1);
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
	
	private Hashtable<LocalDateTime, Double> getEvaporation(String outputFolder) 
																	throws FileNotFoundException
	{
		String outputFile				= outputFolder + "/Aggregated.Values";
		Scanner scanner					= new Scanner(new FileInputStream(new File(outputFile)));
		Hashtable<LocalDateTime, Double> evaporation = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_2);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double value				= Double.valueOf(tokens[8]);
			evaporation.put(dateTime, value);
		}
		scanner.close();
		return evaporation;
	}
	
	private Hashtable<LocalDateTime, double[]> getSoilMoisture(String outputFolder) 
																	throws FileNotFoundException
	{
		String outputFile				= outputFolder + "/Aggregated.Values";
		Scanner scanner					= new Scanner(new FileInputStream(new File(outputFile)));
		Hashtable<LocalDateTime, double[]> soilMoisture = new Hashtable<>();
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_FLOW_2);
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine();
			String[] tokens				= line.split("\\s+");
			LocalDateTime dateTime		= LocalDateTime.parse(tokens[0], formatter);
			double[] sm					= new double[3];
			sm[0]						= Double.valueOf(tokens[30]);
			sm[1]						= Double.valueOf(tokens[31]);
			sm[2]						= Double.valueOf(tokens[32]);
			soilMoisture.put(dateTime, sm);
		}
		scanner.close();
		return soilMoisture;
	}	

	private ArrayList<Double> computeFitness(ArrayList<Double> initialStateArray, 
												Hashtable<LocalDateTime, Double> modeledQ)
	{
		ArrayList<Double> fitness			= new ArrayList<>();
		for (Objective objective : objectives)
		{
			String objID					= objective.getId();
			boolean ok						= true;
			double fitVal					= Double.NaN;
			double[] source					= Utilities.toArray(initialStateArray);
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
					double[] ssArr			= Utilities.toArray(initialStateArray);
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
	
	public void prepareForecast(NonParametric currentState, LocalDateTime dateTime,
								String outputFolder) throws IOException
	{
		// Create mean streamflow file
		meanQFile							= outputFolder + MEAN_Q_FILE;
		PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, false)));
		out.println("Date time\tObserved streamflow (l/s)\tMean streamflow (l/s)\tSt. dev. (l/s)");
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
			State targetState				= config.initialState;
			ArrayList<Soil> soils			= config.soils;
			ArrayList<Vegetation> vegetations = config.vegetations;
			StreamNetwork network			= config.network;
			DHSVMParticle particle			= new DHSVMParticle(this, id, soils, vegetations, 
						network, null, targetArray, null, targetState, null, null);
			ParticleWrapper wrapper			= new ParticleWrapper(particle, null, dateTime);
			wrapper.setWeight(sample.getWeight());
			samples.set(p, wrapper);
		}
	}
	
	public NonParametric forecast(LocalDateTime forecastEnd, String folder, int threadCount, 
						LocalDateTime performanceStart, ArrayList<LocalDateTime> statesToSave,
						boolean removeFolder, long timeLimit) throws IOException
	{
		System.out.println("Starting forecast...");
		this.forecastFolder			= folder;
		this.removeForecastFiles	= removeFolder;
		if (statesToSave == null)
			statesToSave			= new ArrayList<>();
		if (statesToSave.size() == 0)
			statesToSave.add(forecastEnd);
		this.statesToSave			= statesToSave;
		if (!folder.equals(modelsFolder))
			createModelsFolders(folder);
		
		// Prepare distributions of forecasted outputs
		ParticleWrapper particle	= (ParticleWrapper) currentState.getSamples().get(0);
		LocalDateTime forecastStart	= particle.getEnd().plus(modelTimeStep);
		LocalDateTime current		= forecastStart;
		forecastStreamflow			= new Hashtable<>();
		forecastEvaporation			= new Hashtable<>();
		forecastSoilMoistureL1		= new Hashtable<>();
		forecastSoilMoistureL2		= new Hashtable<>();
		forecastSoilMoistureL3		= new Hashtable<>();
		while (!current.isAfter(forecastEnd))
		{
			forecastStreamflow.put(		current, new KernelDensity());
			forecastEvaporation.put(	current, new KernelDensity());
			forecastSoilMoistureL1.put(	current, new KernelDensity());
			forecastSoilMoistureL2.put(	current, new KernelDensity());
			forecastSoilMoistureL3.put(	current, new KernelDensity());
			current					= current.plus(modelTimeStep);
		}
		
		// Verify folders
		Path path					= FileSystems.getDefault().getPath(folder + "/Forecasts");
		if (!Files.exists(path))
			Files.createDirectory(path);
		path						= FileSystems.getDefault().getPath(folder + "/Input");
		if (!Files.exists(path))
			createModelsFolders(folder);
		
		// Populate forecast queue
		int sampleCount				= currentState.getSamples().size();
		forecastEndStates			= new ArrayList<>(sampleCount);
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
					timeUp			= time > timeLimit;
					if (timeUp)
						System.out.println("Forecast time up!");
					wait(1);
				} catch (InterruptedException e) {}
			}
		forecastQueue.clear();
		
		// Compute bandwidth of distribution
		for (LocalDateTime dateTime : forecastStreamflow.keySet())
		{
			forecastStreamflow.get(		dateTime).computeGaussianBandwidth();
			forecastEvaporation.get(	dateTime).computeGaussianBandwidth();
			forecastSoilMoistureL1.get(	dateTime).computeGaussianBandwidth();
			forecastSoilMoistureL2.get(	dateTime).computeGaussianBandwidth();
			forecastSoilMoistureL3.get(	dateTime).computeGaussianBandwidth();
		}
		
		// Write mean streamflow
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(DATE_TIME_FORMAT_REPORT);
		PrintWriter out 	= new PrintWriter(new BufferedWriter(new FileWriter(meanQFile, true)));
		current						= forecastStart;
		while (!current.isAfter(forecastEnd))
		{
			KernelDensity dist		= forecastStreamflow.get(current);
			out.println(formatter.format(current) + "\t" + obsQ.get(current) + "\t" 
									+ dist.getMean() + "\t" + dist.getStDev());
			current					= current.plus(modelTimeStep);
		}
		out.close();
		
		// Compute forecast performance metrics
		if (performanceStart != null)
		{
			// Obtain data
			current					= performanceStart;
			ArrayList<Double> obs	= new ArrayList<>();
			ArrayList<Double> mod	= new ArrayList<>();
			ContSeries density		= new ContSeries(false);
			ContSeries rarity		= new ContSeries(false);
			while (!current.isAfter(forecastEnd))
			{
				KernelDensity dist	= forecastStreamflow.get(current);
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
		
		// Remove file folders
		if (removeFolder)
		{
			try
			{
				FileUtils.deleteDirectory(new File(folder + "/Forecasts"	));
				FileUtils.deleteDirectory(new File(folder + "/Input"		));
				FileUtils.deleteDirectory(new File(folder + "/Met"			));
			} catch (Exception e)
			{
				System.out.println("Error removing folder: " + e.getMessage());
			}
		}
		
		System.out.println("");
		
		// Return final distribution
		if (forecastEndStates.size() > 0)
		{
			try
			{
				synchronized (forecastEndStates)
				{
					return configurator.createDistribution(forecastEndStates);
				}
			} catch (Exception e)
			{
				e.printStackTrace();
				return null;
			}
		}
		else
		{
			System.out.println("No states to create final distribution");
			return null;
		}
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastStreamflow()
	{
		return forecastStreamflow;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastEvaporation()
	{
		return forecastEvaporation;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL1()
	{
		return forecastSoilMoistureL1;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL2()
	{
		return forecastSoilMoistureL2;
	}
	
	public Hashtable<LocalDateTime, KernelDensity> getForecastSoilMoistureL3()
	{
		return forecastSoilMoistureL3;
	}

	@Override
	public void execute(String processID)	// Perform forecast
	{
		while (forecastQueue.size() > 0)
		{
			String runFolder				= null;
			try
			{
				ContMultiSample sample		= forecastQueue.poll();
				LocalDateTime end			= forecastEnd;
				
				// Prepare model
				ParticleWrapper wrapper		= (ParticleWrapper) sample;
				DHSVMParticle particle		= (DHSVMParticle) wrapper.getParticle();
				String particleID			= particle.getId();
				State initialState			= particle.getTargetState();
				LocalDateTime start			= initialState.dateTime;
				runFolder					= forecastFolder + "/Forecasts/" + particleID;
				String stateFolder			= runFolder + "/state";
				String outputFolder			= runFolder + "/output";
				String configFile			= runFolder + "/Configuration.txt";
				Files.createDirectory(FileSystems.getDefault().getPath(runFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(stateFolder	));
				Files.createDirectory(FileSystems.getDefault().getPath(outputFolder	));
				initialState.writeToFiles(stateFolder);
				
				String classFileName				= "stream_class.txt";
				String networkFileName				= "stream_network.txt";
				String mapFileName					= "stream_map.txt";
				String surfaceRoutingFile			= "surface_routing.txt";
				
				ArrayList<Soil> soils				= particle.getSoils();
				ArrayList<Vegetation> vegetations	= particle.getVegetations();
				StreamNetwork network				= particle.getNetwork();
				DHSVMModel model					= new DHSVMModel(defaultModel.optionsSection,
				defaultModel.areaSection, defaultModel.constantsSection, defaultModel.demFile,
				defaultModel.maskFile, defaultModel.flowDirFile, defaultModel.soilTypeFile,
				defaultModel.soilDepthFile, defaultModel.vegTypeFile, classFileName,
				networkFileName, mapFileName, surfaceRoutingFile, defaultModel.d8FlowDirection,
				soils, vegetations, defaultModel.stations);
				
				network.writeToFiles(runFolder, classFileName, networkFileName, mapFileName,
										surfaceRoutingFile);
				
				model.writeConfigFile(configFile, start, end, modelTimeStep, "state", "output", 
										false, statesToSave);
				
				// Run forecast
				ProcessBuilder pb		= new ProcessBuilder(dhsvmExec, configFile);
				pb.directory(new File(runFolder));
				Process process			= pb.start();
				BufferedReader br		= new BufferedReader(new InputStreamReader(
																		process.getInputStream()));
				while (br.readLine() != null)
				{
					synchronized (this)	{try {wait(1);} catch (InterruptedException e) {}}
				}
				
				// Retrieve target outputs
				String line				= particleID + "\t" + sample.getWeight() + "\t";
				Hashtable<LocalDateTime, Double> q		= getHydrograph(	outputFolder);
				Hashtable<LocalDateTime, Double> ev		= getEvaporation(	outputFolder);
				Hashtable<LocalDateTime, double[]> sm	= getSoilMoisture(	outputFolder);
				LocalDateTime current 	= start.plus(modelTimeStep);
				while (!current.isAfter(end))
				{
					line				+= q.get(current) + "\t";
					current				= current.plus(modelTimeStep);
				}
				current					= start.plus(modelTimeStep);
				double weight			= sample.getWeight();
				while (!current.isAfter(forecastEnd))
				{
					Double qi			= q.get(current);
					Double evi			= ev.get(current);
					double[] smc		= sm.get(current);
					if (qi != null)		forecastStreamflow.get(current).addSample(qi, weight);
					if (evi != null)	forecastEvaporation.get(current).addSample(evi, weight);
					if (smc != null)
					{
						forecastSoilMoistureL1.get(current).addSample(smc[0], weight);
						forecastSoilMoistureL2.get(current).addSample(smc[1], weight);
						forecastSoilMoistureL3.get(current).addSample(smc[2], weight);
					}
					current				= current.plus(modelTimeStep);
				}
				
				try
				{
					State finalState		= new State(forecastEnd, runFolder + "/output",
															input.rows, input.cols);
					ModelConfiguration conf	= new ModelConfiguration(soils, vegetations, network,
																		finalState);
					ArrayList<Double> vals	= configurator.toArray(conf);
					forecastEndStates.add(new Sample(weight, vals));
					
					System.out.println("Completed forecast for " + particleID);
				} catch (Exception e)
				{
					System.out.println("Completed partial forecast for "
										+ particleID + ": " + e.getMessage());
				}
				
				forecastQ.put(particleID, line);
			} catch (Exception e)
			{
				e.printStackTrace();
			} finally
			{
				// Delete folder
				if (removeForecastFiles)
					try
					{
						FileUtils.deleteDirectory(new File(runFolder));
					} catch (Exception e)
					{
						// Do nothing
					}
			}
		}
	}

}
