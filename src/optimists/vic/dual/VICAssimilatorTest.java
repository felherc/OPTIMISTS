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
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.TreeSet;

import org.apache.commons.io.FileUtils;

import maestro_mo.ContVar;
import maestro_mo.MAESTRO;
import optimists.OPTIMISTS;
import probDist.ContProbDist;
import probDist.KernelDensity;
import probDist.multiVar.EnsembleGGMLite;
import probDist.multiVar.KD_GGMLite;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.Utilities;
import utilities.geom.Point2D;
import utilities.stat.ContStats;
import vic.Forcing;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICAssimilatorTest
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	public final static String STATES_FOLDER		= "Base states";
	public final static String MODELS_PREP_FOLDER	= "Preparation";
	public final static String DA_FOLDER 			= "Data assimilation";

	public final static String STATE_FILE_DT_FORMAT	= "yyyyMMdd HH-mm";

	public final static String OUT_FILE_STATS		= "Stats.txt";
	public final static String OUT_FILE_PERF		= "Performance.txt";
	public final static String OUT_FILE_Q			= "Q.txt";
	public final static String OUT_FILE_EV			= "Ev.txt";
	public final static String OUT_FILE_SM			= "SM.txt";
	public final static String OUT_FILE_W			= "W.txt";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private VICAssimilator			assimilator;
	
	private Duration				modelTimeStep;
	private boolean					removeDAFiles;
	private boolean					removeDAModelFiles;
	private boolean					removeForecastFiles;
	private String					modelsFolder;
	private int						threadCount;
	
	private String					outputFolder;
	private DateTimeFormatter		formatter;
	
	private int						dimLimit;
	private double					corrThreshold;
	private int						distType;
	private double					scaling;
	private boolean					silverman;
	private GGMLiteCreator			ggmCreator;
	
	private TreeSet<LocalDateTime>	stateList;					
	private String					statesFolder;
	private String					stateFileHeader;
	
	private Hashtable<LocalDateTime, Double> obsQ;
	
	private HashSet<String>			toDelete;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICAssimilatorTest(String parameterFolder, ModelConfigurator configurator,
			ArrayList<String> globalFileParams,	ArrayList<Forcing> cellForcings, 
			Duration modelTimeStep, Hashtable<String, Double> areas, ArrayList<String> outputs,
			boolean defaultParameters, String vicExec, long simMaxTime, int maxAttemptsPerPart,
			int maxForecasts, boolean removeDAFiles, boolean removeDAModelFiles,
			boolean removeForecastFiles, Hashtable<LocalDateTime, Double> obsQ, boolean objQNSE,
			boolean objQMAE, boolean objQMARE, boolean objIndeppdf, boolean objpdf,
			boolean objLogpdf, boolean objMDist, boolean objMForce, boolean objMeanMForce)
					throws IOException
	{
		this.modelTimeStep			= modelTimeStep;
		this.removeDAFiles			= removeDAFiles;
		this.removeDAModelFiles		= removeDAModelFiles;
		this.removeForecastFiles	= removeForecastFiles;
		this.obsQ					= obsQ;
		assimilator					= new VICAssimilator(parameterFolder, configurator,
				globalFileParams, cellForcings, modelTimeStep, areas, outputs, defaultParameters,
				vicExec, simMaxTime, maxAttemptsPerPart, maxForecasts, removeDAModelFiles, obsQ,
				objQNSE, objQMAE, objQMARE, objIndeppdf, objpdf, objLogpdf, objMDist, objMForce,
				objMeanMForce);
		formatter				= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	public void testAssimilator(String testName, int runIndex, String outputFolder,
			String modelsFolder, LocalDateTime forecastStart, LocalDateTime forecastEnd,
			Duration daTimeStep, ArrayList<Duration> leadTimes, NonParametric baseState,
			LocalDateTime baseStateTime, ArrayList<ContVar> variables, MAESTRO maestro,
			int ensembleSize, int candidateCount, int populationSize, int maxEvaluations,
			double samplePercentage, double rootPercentage, int dimLimit, double corrThreshold,
			int distType, double scaling, boolean silverman, GGMLiteCreator ggmCreator,
			boolean weightPerFront, double particleGreed, int threadCount,
			boolean assimilatorReports, boolean maestroReports, String hallOfFameFolder,
			long timeLimit, int rankHistogramBins) throws IOException
	{
		this.outputFolder					= outputFolder;
		this.modelsFolder					= modelsFolder;
		this.threadCount					= threadCount;
		this.dimLimit						= dimLimit;
		this.corrThreshold					= corrThreshold;
		this.distType						= distType;
		this.scaling						= scaling;
		this.silverman						= silverman;
		this.ggmCreator						= ggmCreator;
		
		toDelete							= new HashSet<>();
		
		// Initialize report files
		Collections.sort(leadTimes);
		for (Duration leadTime : leadTimes)
		{
			// Create folder
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;
			Files.createDirectory(FileSystems.getDefault().getPath(ltFolder));
			
			// Create streamflow files
			String mFile					= ltFolder + "/" + OUT_FILE_STATS;
			String qFile					= ltFolder + "/" + OUT_FILE_Q;
			String eFile					= ltFolder + "/" + OUT_FILE_EV;
			String sFile					= ltFolder + "/" + OUT_FILE_SM;
			String wFile					= ltFolder + "/" + OUT_FILE_W;
			PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(mFile, false)));
			out.println("Date-time\tQ_Mean\tQ_stDev\tEv_Mean\tEv_stDev\tSM_Mean\tSM_stDev");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(qFile, false)));
			out.println("Date-time\tStreamflow values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(eFile, false)));
			out.println("Date-time\tEvaporation values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(sFile, false)));
			out.println("Date-time\tSoil moisture values");
			out.close();
			out				= new PrintWriter(new BufferedWriter(new FileWriter(wFile, false)));
			out.println("Date-time\tParticle weights");
			out.close();
		}
		
		// Initialize state folders
		stateList							= new TreeSet<>();
		statesFolder						= outputFolder + "/" + STATES_FOLDER;
		String prepFolder					= modelsFolder + "/" + MODELS_PREP_FOLDER;
		StringBuilder builder				= new StringBuilder("Id\tWeight");
		for (ContVar var : variables)
		{
			builder.append("\t");
			builder.append(var.getName());
		}
		stateFileHeader						= builder.toString();
		Files.createDirectory(FileSystems.getDefault().getPath(prepFolder));
		Files.createDirectory(FileSystems.getDefault().getPath(statesFolder));
		writeStateToFile(baseStateTime, baseState);
		
		// Perform test loop
		Files.createDirectory(FileSystems.getDefault().getPath(outputFolder + "/" + DA_FOLDER));
		LocalDateTime current				= LocalDateTime.MAX;
		LocalDateTime terminate				= LocalDateTime.MIN;
		for (Duration leadTime : leadTimes)
		{
			if (forecastStart.minus(leadTime).isBefore(current))
				current						= forecastStart.minus(leadTime);
			if (forecastEnd.minus(leadTime).isAfter(terminate))
				terminate					= forecastEnd.minus(leadTime);
		}
		current								= current.plus(modelTimeStep);
		while (!current.isAfter(terminate))
		{
			// Prepare assimilation
			LocalDateTime daStart			= current.minus(daTimeStep);
			NonParametric initState			= getState(daStart);
			String strDT					= formatter.format(current);			
			String daRunName				= testName + "_" + strDT;
			String folder					= outputFolder + "/" + DA_FOLDER + "/" + strDT;
			String modelsiFolder			= modelsFolder + "/" + strDT;
			Files.createDirectory(FileSystems.getDefault().getPath(folder));
			Files.createDirectory(FileSystems.getDefault().getPath(modelsiFolder));
			
			// Assimilate and store target state
			assimilator.assimilate(daRunName, runIndex, folder, modelsiFolder, daStart,
					current, daTimeStep, initState, variables, maestro, ensembleSize,
					candidateCount, populationSize, maxEvaluations, samplePercentage,
					rootPercentage, dimLimit, corrThreshold, distType, scaling, silverman,
					ggmCreator, weightPerFront, particleGreed, threadCount, assimilatorReports,
					maestroReports, false, hallOfFameFolder, timeLimit);
			NonParametric targetState		= assimilator.getTargetState();
			writeStateToFile(current, targetState);
			
			// Perform forecast
			LocalDateTime dateTime			= current;
			NonParametric state				= targetState;
			for (Duration leadTime : leadTimes)
			{
				LocalDateTime end			= current.plus(leadTime);
				if (!end.isAfter(forecastStart) || end.isAfter(forecastEnd))
					continue;
				
				strDT						= formatter.format(dateTime);
				assimilator.prepareForecast(state, dateTime, modelsiFolder);
				state						= assimilator.forecast(end, modelsiFolder, false,
												false, threadCount, removeForecastFiles);
				
				KernelDensity streamflow	= assimilator.getFinalStreamflow();
				KernelDensity evaporation	= assimilator.getFinalEvaporation();
				KernelDensity soilMoisture	= assimilator.getFinalSoilMoisture();
				dateTime					= end;
				storeResults(dateTime, leadTime, state, streamflow, evaporation, soilMoisture);
			}
			
			// Delete assimilation model folders
			if (removeDAModelFiles)
			{
				try
				{
					FileUtils.deleteDirectory(new File(modelsiFolder));
				} catch (Exception e)
				{
					toDelete.add(folder);
				}
			}
			
			// Delete folders
			if (removeDAFiles)
			{
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e)
				{
					toDelete.add(folder);
				}
				try
				{
					FileUtils.deleteDirectory(new File(modelsiFolder));
				} catch (Exception e)
				{
					toDelete.add(modelsiFolder);
				}
			}
			
			// Advance to next time step
			current							= current.plus(modelTimeStep);
		}
		
		// Delete remaining folders
		if (removeDAFiles)
		{
			for (String folder : toDelete)
				try
				{
					FileUtils.deleteDirectory(new File(folder));
				} catch (Exception e) {}
		}
		
		// Compute performance
		computeForecastPerformance(leadTimes, forecastStart, forecastEnd, rankHistogramBins);
		
		System.out.println("\nFinished!");
	}
	
	private NonParametric getState(LocalDateTime dateTime) throws IOException
	{
		String strDT				= formatter.format(dateTime);
		
		// Try loading initial state
		try
		{
			NonParametric state		= readStateFromFile(dateTime);
			if (!stateList.contains(dateTime))
				stateList.add(dateTime);
			System.out.println("\nLoaded initial state (" + strDT + ")...");
			return state;
		} catch (Exception e)
		{
			// Continue
		}
		
		System.out.println("\nPreparing initial state for assimilation (" + strDT + ")...");
		LocalDateTime baseDT		= stateList.lower(dateTime);
		NonParametric baseState		= readStateFromFile(baseDT);
		String folder				= modelsFolder + "/" + MODELS_PREP_FOLDER + "/" + strDT;
		Files.createDirectory(FileSystems.getDefault().getPath(folder));
		assimilator.prepareForecast(baseState, baseDT, folder);
		NonParametric state			= assimilator.forecast(dateTime, folder, false, false,
															threadCount, removeDAModelFiles);
		writeStateToFile(dateTime, state);
		stateList.add(dateTime);
		return state;
	}

	private void writeStateToFile(LocalDateTime dateTime, NonParametric stateDistribution)
					throws IOException
	{
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(fileName, false)));
		out.println(stateFileHeader);
		int index					= 1;
		for (ContMultiSample particle : stateDistribution.getSamples())
		{
			StringBuilder builder	= new StringBuilder(index + "\t" + particle.getWeight());
			for (Double value : particle.getValues())
			{
				builder.append("\t");
				builder.append(value);
			}
			out.println(builder.toString());
		}
		out.close();
		stateList.add(dateTime);
	}
	
	private NonParametric readStateFromFile(LocalDateTime dateTime)
	{
		DateTimeFormatter formatter	= DateTimeFormatter.ofPattern(STATE_FILE_DT_FORMAT);
		String fileName				= statesFolder + "/" + formatter.format(dateTime) + ".txt";
		Scanner scanner				= null;
		try
		{
			// Read samples from file
			scanner					= new Scanner(new FileInputStream(new File(fileName)));
			scanner.nextLine();
			ArrayList<ContMultiSample> samples = new ArrayList<>();
			while (scanner.hasNextLine())
			{
				String[] tokens		= scanner.nextLine().split("\t");
				double weight		= Double.valueOf(tokens[1]);
				int valueCount		= tokens.length - 2;
				ArrayList<Double> values = new ArrayList<>(valueCount);
				for (int v = 2; v < valueCount + 2; v++)
					values.add(Double.valueOf(tokens[v]));
				samples.add(new Sample(weight, values));
			}
			
			// Create distribution
			NonParametric dist		= null;
			if (distType == OPTIMISTS.TYPE_D_KERNEL || distType == OPTIMISTS.TYPE_F_KERNEL)
			{
				dist				= new MultiVarKernelDensity(true);
				dist.setSamples(samples);
				if (distType == OPTIMISTS.TYPE_D_KERNEL)
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
			else if (distType == OPTIMISTS.TYPE_GGM_LITE)
			{
				dist				= new EnsembleGGMLite(true);
				dist.setSamples(samples);
			}
			else if (distType == OPTIMISTS.TYPE_KD_GGM_LITE)
			{
				dist				= new KD_GGMLite(true);
				dist.setSamples(samples);
				if (ggmCreator == null)
					((KD_GGMLite)dist).computeGaussianBW(scaling);
				else
					((KD_GGMLite)dist).computeGaussianBW(scaling, ggmCreator);
			}
			return dist;
		} catch (Exception e)
		{
			throw new RuntimeException(e.getMessage());
		} finally
		{
			scanner.close();
		}
	}

	private void storeResults(LocalDateTime dateTime, Duration leadTime,
			NonParametric state, KernelDensity streamflow, KernelDensity evaporation,
			KernelDensity soilMoisture) throws IOException
	{
		// Get strings
		String ltStr		= leadTime.toString();
		String ltFolder		= outputFolder + "/Lead time = " + ltStr;
		String dtStr		= formatter.format(dateTime);
		
		String mFile		= ltFolder + "/" + OUT_FILE_STATS;
		String qFile		= ltFolder + "/" + OUT_FILE_Q;
		String eFile		= ltFolder + "/" + OUT_FILE_EV;
		String sFile		= ltFolder + "/" + OUT_FILE_SM;
		String wFile		= ltFolder + "/" + OUT_FILE_W;
		
		// Write statistics
		Path mPath			= Paths.get(mFile);
		String line			= System.lineSeparator() + dtStr + "\t"
									+ streamflow.getMean() + "\t" + streamflow.getStDev() + "\t"
									+ evaporation.getMean() + "\t" + evaporation.getStDev() + "\t"
									+ soilMoisture.getMean() + "\t" + soilMoisture.getStDev();
		BufferedWriter writer = Files.newBufferedWriter(mPath, StandardOpenOption.APPEND);
		writer.write(line);
		
		// Write streamflow and weight values
		PrintWriter outQ	= new PrintWriter(new BufferedWriter(new FileWriter(qFile, true)));
		PrintWriter outW	= new PrintWriter(new BufferedWriter(new FileWriter(wFile, true)));
		String lineQ		= dtStr;
		String lineW		= dtStr;
		ArrayList<Point2D> samples = streamflow.getSamplesWeights();
		for (Point2D sample : samples)
		{
			lineQ			+= "\t" + sample.x;
			lineW			+= "\t" + sample.y;
		}
		outQ.println(lineQ);
		outW.println(lineW);
		outQ.close();
		outW.close();
		
		// Write evaporation values
		PrintWriter outE	= new PrintWriter(new BufferedWriter(new FileWriter(eFile, true)));
		String lineE		= dtStr;
		samples				= evaporation.getSamplesWeights();
		for (Point2D sample : samples)
			lineE			+= "\t" + sample.x;
		outE.println(lineE);
		outE.close();
		
		// Write soil moisture values
		PrintWriter outS	= new PrintWriter(new BufferedWriter(new FileWriter(sFile, true)));
		String lineS		= dtStr;
		samples				= soilMoisture.getSamplesWeights();
		for (Point2D sample : samples)
			lineS			+= "\t" + sample.x;
		outS.println(lineS);
		outS.close();
	}

	private void computeForecastPerformance(ArrayList<Duration> leadTimes,
			LocalDateTime forecastStart, LocalDateTime forecastEnd, int bins) throws IOException
	{
		System.out.println("\nForecast performance:");
		System.out.println("Lead time\tNSE_l2\tNSE_l1\tMARE\tCRPS\tDensity\tHistogram");
		for (Duration leadTime : leadTimes)
		{
			String ltStr					= leadTime.toString();
			String ltFolder					= outputFolder + "/Lead time = " + ltStr;
			String valuesFile				= ltFolder + "/" + OUT_FILE_Q;
			String weightsFile				= ltFolder + "/" + OUT_FILE_W;
			ArrayList<ContProbDist> series	= new ArrayList<>();
			Scanner scannerV				= new Scanner(new FileInputStream(new File(valuesFile)));
			Scanner scannerW				= new Scanner(new FileInputStream(new File(weightsFile)));
			scannerV.nextLine();	// Header
			scannerW.nextLine();	// Header
			
			// Load distribution time series from files
			ArrayList<Double> meanQ			= new ArrayList<>();
			ArrayList<Double> obs			= new ArrayList<>();
			int[] rankHistogram				= new int[bins]; 
			while (scannerV.hasNextLine())
			{
				String[] tokensV			= scannerV.nextLine().split("\t");
				String[] tokensW			= scannerW.nextLine().split("\t");
				
				LocalDateTime dateTime		= LocalDateTime.parse(tokensV[0], formatter);
				Double obsi					= obsQ.get(dateTime);
				if (obsi != null 	&& dateTime.isAfter(	forecastStart	)
									&& !dateTime.isAfter(	forecastEnd		))
				{
					KernelDensity dist		= new KernelDensity();
					dist.setWeighted(true);
					for (int v = 1; v < Math.min(tokensV.length, tokensW.length); v++)
					{
						double value		= Double.valueOf(tokensV[v]);
						double weight		= Double.valueOf(tokensW[v]);
						dist.addSample(value, weight);
					}
					if (dist.getSamples().size() == 0)
					{
						scannerV.close();
						scannerW.close();
						throw new RuntimeException("No forecasted values at time " + tokensV[0]);
					}
					dist.computeGaussianBandwidth();
					double quantile			= dist.getCDF(obsi);
					int bin					= Math.min(bins - 1, (int)Math.floor((bins*quantile)));
					rankHistogram[bin]++;
					series.add(dist);
					meanQ.add(dist.getMean());
					obs.add(obsi);
				}
			}
			scannerV.close();
			scannerW.close();
			
			// Compute density and CRPS
			ContStats density				= new ContStats(false);
			ContStats crps					= new ContStats(false);
			for (int t = 0; t < obs.size(); t++)
			{
				double obsi					= obs.get(t);
				ContProbDist dist			= series.get(t);
				density.addValue(dist.getpdf(obsi));
				KernelDensity dist2			= (KernelDensity)dist;
				crps.addValue(dist2.computeEnsembleCRPS(obsi));
			}
			
			// Compute rank histogram heterogeneity
			
			double average					= ((double)series.size())/bins;
			ContStats deviation				= new ContStats(false);
			for (int b = 0; b < bins; b++)
				deviation.addValue(Math.abs(rankHistogram[b] - average)/average);
			double rankHistHetero			= deviation.getMean();

			// Compute performance
			double[] obsArr					= Utilities.toArray(obs);
			double[] modArr					= Utilities.toArray(meanQ);
			double nse_l2					= Utilities.computeNashSutcliffe(obsArr, modArr);
			double nse_l1					= Utilities.computeNashSutcliffe(obsArr, modArr, 1.0);
			double mare						= Utilities.computeMeanAbsRelativeError(obsArr, modArr);
			double meanDensity				= density.getMean();
			double meanCRPS					= crps.getMean();

			// Store performance file
			String performanceFile	= ltFolder + "/" + OUT_FILE_PERF;
			PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(performanceFile)));
			out.println("Values\t"						+ obsArr.length		);
			out.println("NSE_l2\t"						+ nse_l2			);
			out.println("NSE_l1\t"						+ nse_l1			);
			out.println("MARE\t"						+ mare				);
			out.println("Mean CRPS\t"					+ meanCRPS			);
			out.println("Mean density\t"				+ meanDensity		);
			out.println("Rank history heterogeneity\t"	+ rankHistHetero	);
			out.println("");
			out.println("Rank histogram:");
			out.println("");
			out.println("Bin\tCount");
			for (int b = 0; b < bins; b++)
				out.println((b + 1) + "\t" + rankHistogram[b]);
			out.close();
			
			// Print to console
			System.out.println(ltStr + "\t" + nse_l2 + "\t" + nse_l1 + "\t" + mare + "\t"
								+ meanCRPS + "\t" + meanDensity + "\t" + rankHistHetero);
		}
	}
	
}
