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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;

import maestro_mo.ContVar;
import maestro_mo.MAESTRO;
import maestro_mo.gen.GA;
import maestro_mo.gen.MetroACO;
import optimists.OPTIMISTS;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumElement;
import vic.routing.MuskingumNetwork;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class BlueRiverTest
{

	public static void main(String[] args) throws IOException
	{
		// General
		String testName					= "Blue River OPTIMISTS test";
		int runIndex					= 1;
		String outputFolder				= "data/Test " + runIndex;
		String modelsFolder				= outputFolder + "/Models";
		String inputDataFolder			= "data/Blue River";
		String hallOfFameFolder			= "";
		boolean assimilatorReports		= false;
		boolean maestroReports			= false;
		boolean defaultParameters		= false;
		boolean removeDAFiles			= false;
		boolean removeDAModelFiles		= true;
		boolean removeForecastFiles		= true;
		
		// Determine times
		LocalDateTime forecastStart		= LocalDateTime.of(1996, 12, 31, 0, 0);
		LocalDateTime forecastEnd		= LocalDateTime.of(1997,  6,  1, 0, 0);
		Duration modelTimeStep			= Duration.ofDays(1);
		Duration daTimeStep				= Duration.ofDays(7);
		ArrayList<Duration> leadTimes	= new ArrayList<>();
		leadTimes.add(Duration.ofDays(1));
		leadTimes.add(Duration.ofDays(3));
		leadTimes.add(Duration.ofDays(6));
		leadTimes.add(Duration.ofDays(12));
		LocalDateTime baseStateTime		= LocalDateTime.of(1996, 12, 1, 0, 0);
		String initStateFile			= inputDataFolder + "/States/19961201 00-00.txt";
		int timeLimit					= 5*24*60*60*1000;
		int rankHistogramBins			= 10;
		
		// Determine objectives
		boolean objQNSE					= true;
		boolean objQMAE					= false;
		boolean objQMARE				= true;
		boolean objIndeppdf				= true;
		boolean objpdf					= false;
		boolean objLogpdf				= false;
		boolean objMDist				= false;
		boolean objMForce				= false;
		boolean objMeanMForce			= false;
		
		// Determine OPTIMISTS parameters
		int ensembleSize				= 30;
		int candidateCount				= 30;
		int populationSize				= 15;
		int maxEvaluations				= 90;
		double samplePercentage			= 1.0;
		double rootPercentage			= 0.8;
		double randomSolutionRatio		= 0.0;
		int dimLimit					= Integer.MAX_VALUE;
		double corrThreshold			= 0.0;
		boolean weightPerFront			= true;
		double particleGreed			= 0.4;
		int distType					= OPTIMISTS.TYPE_D_KERNEL;
		double scaling					= 0.1;
		boolean silverman				= true;
		GGMLiteCreator ggmCreator		= null;
		int threadCount					= 1;
		
		// Simulation options
		String vicExec					= "data/VIC/vicNl.exe";
		long simMaxTime					= 6*1000;
		int maxEvaluationsPerParticle	= 3;
		int maxForecasts				= 3;
		
		// Read global parameters file and forcings
		String parameterFolder			= inputDataFolder + "/Parameters";
		String soilFile					= parameterFolder + "/soil.txt";
		ArrayList<Soil> soils			= Soil.readFromFile(soilFile, 3);
		String paramFile				= inputDataFolder + "/vic_global_file_val";
		ArrayList<String> globalFileParams = new ArrayList<>();
		Scanner scanner					= new Scanner(new FileInputStream(new File(paramFile)));
		while (scanner.hasNextLine())
			globalFileParams.add(scanner.nextLine());
		scanner.close();
		ArrayList<Forcing> cellForcings	= Forcing.loadFromFiles(inputDataFolder + "/Forcing");
		
		// Read network file
		String routingFile				= parameterFolder + "/routing.txt";
		Hashtable<String, Double> areas = new Hashtable<>();
		HashSet<String> inNetwork		= new HashSet<>();
		ArrayList<String> outputs		= new ArrayList<>();
		Hashtable<String, Double> directFractions = new Hashtable<>();
		MuskingumNetwork network		= new MuskingumNetwork();
		scanner							= new Scanner(new FileInputStream(new File(routingFile)));
		while (scanner.hasNextLine())
		{
			String line					= scanner.nextLine().trim();
			String[] tokens				= line.split("\t");
			String id					= tokens[0];
			double k					= Double.valueOf(tokens[1]);
				double x				= Double.valueOf(tokens[2]);
			String downstream			= null;
			if (tokens.length > 3)
			{
				downstream				= tokens[3];
				if (tokens.length > 4)
				{
					double area			= Double.valueOf(tokens[4]);
					areas.put(id, area);
					if (tokens.length > 5)
					{
						directFractions.put(id, Double.valueOf(tokens[5]));
						if (tokens.length > 6)
							if (Boolean.valueOf(tokens[6]))
								outputs.add(id);
					}
				if (!inNetwork.contains(downstream))
					downstream			= null;
				}
			}
			network.addElement(id, new MuskingumElement(k, x, 0.0, 0.0), downstream);
			inNetwork.add(id);
		}
		scanner.close();
		
		// Setup initial distribution
		String exStateFolder			= inputDataFolder + "/states/scenario 1/";
		State exampleState = new State(exStateFolder + "01", exStateFolder + "01_rout.txt");
		ArrayList<ContVar> variables	= getVariables(soils, exampleState);		
		BlueRiver configurator			= new BlueRiver(soils, network, exampleState,
												directFractions, dimLimit, corrThreshold,
												distType, scaling, silverman, ggmCreator);
		NonParametric baseState 		= createInitialState(initStateFile, variables,
												ensembleSize, distType, scaling, silverman,
												dimLimit, ggmCreator, corrThreshold, configurator);
		
		// Initialize optimizer
		MAESTRO maestro					= new MAESTRO("", 0, null, null, false, false);
		maestro.setRandomSolutionRatio(randomSolutionRatio);
		if (samplePercentage < 1.0)
		{
			GA ga1						= new GA();
			ga1.setRandomMutation(0.1);
			GA ga2						= new GA();
			ga2.setRandomMutation(0.0);
			MetroACO metroACO1			= new MetroACO();
			metroACO1.setUniformProb(0.01);
			MetroACO metroACO2			= new MetroACO();
			metroACO2.setUniformProb(0.0);
			maestro.addGenerator(ga1);
			maestro.addGenerator(ga2);
			maestro.addGenerator(metroACO1);
			maestro.addGenerator(metroACO2);
		}
		
		// Load observed flow
		LocalDateTime qObsStart			= LocalDateTime.of(1996, 10, 16, 0, 0);
		Hashtable<LocalDateTime, Double> obsQ = loadObsHydrograph(inputDataFolder + "/obsQ.txt",
																	qObsStart, modelTimeStep);
		
		VICAssimilatorTest test			= new VICAssimilatorTest(parameterFolder, configurator,
				globalFileParams, cellForcings, modelTimeStep, areas, outputs, defaultParameters,
				vicExec, simMaxTime, maxEvaluationsPerParticle, maxForecasts, removeDAFiles,
				removeDAModelFiles, removeForecastFiles, obsQ, objQNSE, objQMAE, objQMARE,
				objIndeppdf, objpdf, objLogpdf, objMDist, objMForce, objMeanMForce);
		test.testAssimilator(testName, runIndex, outputFolder, modelsFolder, forecastStart,
				forecastEnd, daTimeStep, leadTimes, baseState, baseStateTime, variables,
				maestro, ensembleSize, candidateCount, populationSize, maxEvaluations,
				samplePercentage, rootPercentage, dimLimit, corrThreshold, distType, scaling,
				silverman, ggmCreator, weightPerFront, particleGreed, threadCount,
				assimilatorReports, maestroReports, hallOfFameFolder, timeLimit,
				rankHistogramBins);
	}
	
	private static ArrayList<ContVar> getVariables(ArrayList<Soil> soils, State exampleState)
	{
		ArrayList<ContVar> variables	= new ArrayList<>();
		
		// Add global variables
		variables.add(new ContVar("ws",							1E-5,	    1.0));
		variables.add(new ContVar("c",							 0.1,	    5.0));
		variables.add(new ContVar("ksatMult",					-0.9,	    0.9));
		variables.add(new ContVar("bubble",						 0.1,	  100.0));
		variables.add(new ContVar("bulkMult",					-0.9,	    0.9));
		variables.add(new ContVar("soil/bulk",					 1.01,	    3.0));
		variables.add(new ContVar("wrcMult",					-0.9,	    0.9));
		variables.add(new ContVar("rough",						1E-3,	    0.1));
		variables.add(new ContVar("snow_rough",					1E-3,	    0.1));
		variables.add(new ContVar("resid_moist",				1E-3,	    0.1));
		variables.add(new ContVar("cellX",						1E-3,	    0.5));
		variables.add(new ContVar("channelX",					1E-3,	    0.5));
		variables.add(new ContVar("channelKMult",				 0.01,	    5.0));
		
		// Add cell variables
		for (int s = 1; s <= soils.size(); s++)
		{
			variables.add(new ContVar(s + "_infilt",			 0.0,	    5.0));
			variables.add(new ContVar(s + "_ds",				1E-5,	    1.0));
			variables.add(new ContVar(s + "_dsmax",				1E-5,	   50.0));
			variables.add(new ContVar(s + "_expt",				 0.0,	    5.0));
			variables.add(new ContVar(s + "_exptMult",			-0.9,	    0.9));
			variables.add(new ContVar(s + "_ksat",				 0.0,	10000.0));
			variables.add(new ContVar(s + "_bulk_density",	  1000.0,	 2000.0));
			variables.add(new ContVar(s + "_wrc_fract",			 0.01,	    0.6));
			variables.add(new ContVar(s + "_wpwp/wrc",			 0.01,	    1.0));
			variables.add(new ContVar(s + "_routingK",		 	 0.01,      5.0));
			variables.add(new ContVar(s + "_directFraction",    0.0,	    1.0));
		}
		
		// Add state variables and return
		variables.addAll(exampleState.getVariables());
		return variables;
	}
	
	private static NonParametric createInitialState(String initStateFile,
				ArrayList<ContVar> variables, int sampleCount, int distType, double scaling, 
				boolean silverman, int dimLimit, GGMLiteCreator ggmCreator, double corrThreshold,
				ModelConfigurator configurator) throws FileNotFoundException
	{
		Scanner scanner				= new Scanner(new FileInputStream(new File(initStateFile)));
		String line					= scanner.nextLine();
		
		// Verify variables match
		String[] tokens				= line.split("\t");
		for (int i = 0; i < variables.size(); i++)
		{
			String varName			= variables.get(i).getName();
			String onFile			= tokens[i + 2];
			if (varName.compareTo(onFile) != 0)
			{
				scanner.close();
				throw new IllegalArgumentException("Variable name mismatch: " 
													+ varName + ", " + onFile);
			}
		}
		
		// Import initial state
		ArrayList<ContMultiSample> samples = new ArrayList<>();
		while (scanner.hasNextLine())
		{
			line					= scanner.nextLine();
			tokens					= line.split("\t");
			double weight			= Double.valueOf(tokens[1]);
			ArrayList<Double> values = new ArrayList<>();
			for (int t = 2; t < tokens.length; t++)
				values.add(Double.valueOf(tokens[t]));
			samples.add(new Sample(weight, values));
		}
		scanner.close();
		
		// Create distribution
		Collections.shuffle(samples);
		int origSize						= samples.size();
		for (int s = sampleCount; s < origSize; s++)
			samples.remove(samples.size() - 1);
		return configurator.createDistribution(samples);
	}

	private static Hashtable<LocalDateTime, Double> loadObsHydrograph(String fileRoute, 
			LocalDateTime start, Duration timeStep) throws FileNotFoundException
	{
		Hashtable<LocalDateTime, Double> qObs = new Hashtable<>();
		Scanner scanner					= new Scanner(new FileInputStream(new File(fileRoute)));
		LocalDateTime dateTime			= start;
		while (scanner.hasNextLine())
		{
			Double value				= Double.valueOf(scanner.nextLine());
			qObs.put(dateTime, value);
			dateTime					= dateTime.plus(timeStep);
		}
		scanner.close();
		return qObs;
	}
	
}
