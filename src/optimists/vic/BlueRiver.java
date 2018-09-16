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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;

import maestro_mo.MAESTRO;
import maestro_mo.gen.GA;
import maestro_mo.gen.MetroACO;
import optimists.OPTIMISTS;
import probDist.multiVar.tools.GGMLiteCreator;
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
public class BlueRiver
{

	public static void main(String[] args) throws IOException
	{
		// General configuration
		int scenario					= 1;
		int runIndex					= 1;
		String problemName				= "Blue River data assimilation, scenario " + scenario;
		String outputFolder				= "data/Tests/Scenario " + scenario + "/0" + runIndex;
		String modelsFolder				= outputFolder + "/Models";
		String forecastFolder			= outputFolder;
		String inputDataFolder			= "data/Blue River/";
		boolean assimilatorReports		= true;
		boolean maestroReports			= false;
		String hallOfFameFolder			= "";
		boolean computePerformance		= true;
		
		// Determine objectives
		boolean objQNSE					= false;
		boolean objQMAE					= true;
		boolean objQMARE				= false;
		boolean objIndeppdf				= false;
		boolean objpdf					= false;
		boolean objLogpdf				= false;
		boolean objMDist				= false;
		boolean objMForce				= false;
		boolean objMeanMForce			= true;
		
		// Determine OPTIMISTS parameters
		int ensembleSize				= 30;
		int candidateCount				= 100;
		int populationSize				= 25;
		double samplePercentage			= 1.0;
		double rootPercentage			= 0.95;
		double particleGreed			= 0.75;
		int distributionType			= OPTIMISTS.TYPE_F_KERNEL;
		double scaling					= Double.NaN;
		boolean silverman				= true;
		GGMLiteCreator ggmCreator		= null;
		int threadCount					= 6;
		
		// Determine times
		Duration modelTimeStep			= Duration.ofDays(1);
		Duration daTimeStep				= Duration.ofDays(5);
		int timeLimit					= 12*60*60*1000;
		//int timeLimit					= 120*60*1000;
		LocalDateTime start				= null;
		LocalDateTime end				= null;
		LocalDateTime forecastEnd		= null;
		if (scenario == 1)
		{
			start						= LocalDateTime.of(1996, 10, 15, 0, 0);
			end							= LocalDateTime.of(1996, 10, 29, 0, 0);
			forecastEnd					= LocalDateTime.of(1996, 11, 12, 0, 0);
		}
		else if (scenario == 2)
		{
			start						= LocalDateTime.of(1997,  1, 15, 0, 0);
			end							= LocalDateTime.of(1997,  1, 29, 0, 0);
			forecastEnd					= LocalDateTime.of(1997,  2, 12, 0, 0);
		}
		else if (scenario == 3)
		{
			start						= LocalDateTime.of(1997,  6,  1, 0, 0);
			end							= LocalDateTime.of(1997,  6, 15, 0, 0);
			forecastEnd					= LocalDateTime.of(1997,  6, 29, 0, 0);
		}
		
		// Simulation options
		String vicExec					= "data/VIC/vicNl.exe";
		long simMaxTime					= 45*1000;
		boolean removeFiles				= true;
		
		// Read global parameters file and forcings
		String parameterFolder			= inputDataFolder + "Parameters";
		String soilFile					= parameterFolder + "/soil.txt";
		ArrayList<Soil> soils			= Soil.readFromFile(soilFile, 3);
		String paramFile				= inputDataFolder + "vic_global_file_val";
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
		
		// Initial states
		ArrayList<State> initialStates	= new ArrayList<>();
		String initStateFolder			= inputDataFolder + "/States/scenario " + scenario + "/";
		initialStates.add(new State(initStateFolder + "01", initStateFolder + "01_rout.txt"));
		initialStates.add(new State(initStateFolder + "02", initStateFolder + "02_rout.txt"));
		initialStates.add(new State(initStateFolder + "03", initStateFolder + "03_rout.txt"));
		initialStates.add(new State(initStateFolder + "04", initStateFolder + "04_rout.txt"));
		initialStates.add(new State(initStateFolder + "05", initStateFolder + "05_rout.txt"));
		initialStates.add(new State(initStateFolder + "06", initStateFolder + "06_rout.txt"));
		initialStates.add(new State(initStateFolder + "07", initStateFolder + "07_rout.txt"));
		initialStates.add(new State(initStateFolder + "08", initStateFolder + "08_rout.txt"));
		
		// Initialize optimizer
		MAESTRO maestro					= new MAESTRO("", 0, null, null, false, true);
		GA ga1							= new GA();
		ga1.setRandomMutation(0.1);
		GA ga2							= new GA();
		ga2.setRandomMutation(0.0);
		MetroACO metroACO1				= new MetroACO();
		metroACO1.setUniformProb(0.01);
		MetroACO metroACO2				= new MetroACO();
		metroACO2.setUniformProb(0.0);
		maestro.addGenerator(ga1);
		maestro.addGenerator(ga2);
		maestro.addGenerator(metroACO1);
		maestro.addGenerator(metroACO2);
		
		// Load observed flow
		LocalDateTime qObsStart			= start.plus(modelTimeStep);
		Hashtable<LocalDateTime, Double> obsQ = loadObsHydrograph(inputDataFolder + "obsQ_" 
													+ scenario + ".txt", qObsStart, modelTimeStep);
		
		// Perform assimilation
		VICAssimilator assimilator		= new VICAssimilator(parameterFolder, soils, 
				globalFileParams, cellForcings, modelTimeStep, network, areas, outputs, 
				directFractions, vicExec, simMaxTime, removeFiles, obsQ, objQNSE, objQMAE, 
				objQMARE, objIndeppdf, objpdf, objLogpdf, objMDist, objMForce, objMeanMForce);
		assimilator.assimilate(problemName, runIndex, outputFolder, modelsFolder, start, end, 
				daTimeStep, maestro, initialStates, ensembleSize, candidateCount, populationSize, 
				samplePercentage, rootPercentage, distributionType, scaling, silverman, ggmCreator,
				particleGreed, threadCount, assimilatorReports, maestroReports, hallOfFameFolder,
				timeLimit);
		assimilator.forecast(forecastEnd, forecastFolder, threadCount, computePerformance);
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
