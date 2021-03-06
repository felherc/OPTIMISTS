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

package optimists.vic.dual;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;

import maestro_mo.ContVar;
import optimists.OPTIMISTS;
import probDist.KernelDensity;
import probDist.multiVar.KD_GGMLite;
import probDist.multiVar.MultiVarKernelDensity;
import probDist.multiVar.NonParametric;
import probDist.multiVar.tools.ContMultiSample;
import probDist.multiVar.tools.GGMLiteCreator;
import probDist.multiVar.tools.Sample;
import utilities.geom.Point2D;
import vic.Forcing;
import vic.Soil;
import vic.routing.MuskingumElement;
import vic.routing.MuskingumNetwork;
import vic.routing.State;

/**
 * Citation: Hernández, F. and Liang, X.: Hybridizing Bayesian and variational data assimilation
 * for high-resolution hydrologic forecasting, Hydrol. Earth Syst. Sci., 22, 5759-5779,
 * https://doi.org/10.5194/hess-22-5759-2018, 2018.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class BlueRiverSimulation implements ModelConfigurator
{

	public static void main(String[] args) throws IOException
	{
		// General configuration
		String outputFolder				= "data/Tests/Output";
		String modelsFolder				= outputFolder + "/Models";
		String forecastFolder			= outputFolder;
		String inputDataFolder			= "data/Blue River/";
		String exampleStateFile			= inputDataFolder + "States/scenario 1/01";
		String exampleRoutingFile		= inputDataFolder + "States/scenario 1/01_rout.txt";
		String initStateFile			= inputDataFolder + "States/19970131 00-00.txt";
		boolean defaultParameters		= true;
		
		// Determine OPTIMISTS parameters
		int ensembleSize				= 100;
		int dimLimit					= Integer.MAX_VALUE;
		double corrThreshold			= 0.0;
		int distributionType			= OPTIMISTS.TYPE_F_KERNEL;
		double scaling					= 0.01;
		boolean silverman				= true;
		GGMLiteCreator ggmCreator		= null;
		int threadCount					= 1;
		
		// Determine times
		Duration modelTimeStep			= Duration.ofDays(1);
		LocalDateTime start				= LocalDateTime.of(1997, 1, 31, 0, 0);
		LocalDateTime end				= LocalDateTime.of(1997, 5,  1, 0, 0);
		
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
		
		// Simulation options
		String vicExec					= "data/VIC/vicNl.exe";
		long simMaxTime					= 30*1000;
		int maxEvalsPerParticle			= 2;
		int maxForecasts				= 3;
		boolean removeFiles				= false;
		
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
		State exampleState				= new State(exampleStateFile, exampleRoutingFile);
		ArrayList<ContVar> variables	= getVariables(soils, exampleState);
		BlueRiverSimulation configurator = new BlueRiverSimulation(soils, network, exampleState,
												directFractions, dimLimit, corrThreshold,
												distributionType, scaling, silverman, ggmCreator);
		ArrayList<ContMultiSample> initialState = createInitialState(initStateFile, variables,
											ensembleSize, distributionType, scaling, silverman,
											dimLimit, ggmCreator, corrThreshold, configurator);
		
		// Perform forecast
		VICAssimilator assimilator		= new VICAssimilator(parameterFolder, configurator, 
				globalFileParams, cellForcings, modelTimeStep, areas, outputs, defaultParameters, 
				vicExec, simMaxTime, maxEvalsPerParticle, maxForecasts, removeFiles, null, false,
				false, false, false, false, false, false, false, false);
		assimilator.prepareForecast(initialState, start, modelsFolder);
		ArrayList<ContMultiSample> finalState = assimilator.forecast(end, forecastFolder, true,
														false, threadCount, false);
		Hashtable<LocalDateTime, KernelDensity> qts = assimilator.getDistForecastQ();
		Hashtable<LocalDateTime, KernelDensity> sts = assimilator.getDistForecastSM();
		
		// Write final state file
		StringBuilder builder			= new StringBuilder("Id\tWeight");
		for (ContVar var : variables)
		{
			builder.append("\t");
			builder.append(var.getName());
		}
		DateTimeFormatter formatter		= DateTimeFormatter.ofPattern("yyyyMMdd HH-mm");
		String stateFile				= outputFolder + "/" + formatter.format(end) + ".txt";
		PrintWriter out	= new PrintWriter(new BufferedWriter(new FileWriter(stateFile, false)));
		out.println(builder.toString());
		int index					= 1;
		for (ContMultiSample particle : finalState)
		{
			builder						= new StringBuilder(index + "\t" + particle.getWeight());
			for (Double value : particle.getValues())
			{
				builder.append("\t");
				builder.append(value);
			}
			out.println(builder.toString());
		}
		out.close();
		
		// Write output files
		String qFile					= outputFolder + "/Q.txt";
		String wFile					= outputFolder + "/W.txt";
		String sFile					= outputFolder + "/SM.txt";
		PrintWriter outQ = new PrintWriter(new BufferedWriter(new FileWriter(qFile)));
		PrintWriter outW = new PrintWriter(new BufferedWriter(new FileWriter(wFile)));
		PrintWriter outS = new PrintWriter(new BufferedWriter(new FileWriter(sFile)));
		outQ.println("Date time\tStreamflow values");
		outW.println("Date time\tWeights");
		outS.println("Date time\tSoil moisture");
		LocalDateTime current			= start.plus(modelTimeStep);
		while (!current.isAfter(end))
		{
			String lineQ				= current.toString();
			String lineW				= current.toString();
			String lineS				= current.toString();
			KernelDensity distQ			= qts.get(current);
			for (Point2D point : distQ.getSamplesWeights())
			{
				lineQ					+= "\t" + point.x;
				lineW					+= "\t" + point.y;
			}
			KernelDensity distS			= sts.get(current);
			for (Point2D point : distS.getSamplesWeights())
				lineS					+= "\t" + point.x;
			outQ.println(lineQ);
			outW.println(lineW);
			outS.println(lineS);
			current						= current.plus(modelTimeStep);
		}
		outQ.close();
		outW.close();
		outS.close();
		
		System.out.println("\nFinished!");
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
	
	private static ArrayList<ContMultiSample> createInitialState(String initStateFile,
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
		return samples;
	}
	
	private ArrayList<Soil>				origSoils;
	private MuskingumNetwork			origNetwork;
	private State						origState;
	private Hashtable<String, Double>	origDirectFractions;
	
	private int							dimLimit;
	private double						corrThreshold;
	private int							distType;
	private double						scaling;
	private boolean						silverman;
	private GGMLiteCreator				ggmCreator;
	
	public BlueRiverSimulation(ArrayList<Soil> origSoils, MuskingumNetwork origNetwork, State origState,
						Hashtable<String, Double> origDirectFractions, int dimLimit,
						double corrThreshold, int distType, double scaling, boolean silverman,
						GGMLiteCreator ggmCreator)
	{
		this.origSoils					= origSoils;
		this.origNetwork				= origNetwork;
		this.origState					= origState;
		this.origDirectFractions		= origDirectFractions;
		this.dimLimit					= dimLimit;
		this.corrThreshold				= corrThreshold;
		this.distType					= distType;
		this.scaling					= scaling;
		this.silverman					= silverman;
		this.ggmCreator					= ggmCreator;
	}

	public ModelConfiguration configure(ArrayList<Double> values, LocalDateTime dateTime, 
										boolean defaultParameters)
	{
		// Clone original values
		ArrayList<Soil> soils			= Soil.cloneList(origSoils);
		MuskingumNetwork network		= origNetwork.clone();
		Hashtable<String, Double> directFractions = origDirectFractions;
		
		// Read state
		int parameterCount				= 13 + soils.size()*11;
		ArrayList<Double> stateArray	= new ArrayList<>(values.size() - parameterCount);
		for (int s = parameterCount; s < values.size(); s++)
			stateArray.add(values.get(s));
		State state						= new State(stateArray, dateTime, origState);
		
		if (defaultParameters)
			return new ModelConfiguration(state, soils, network, directFractions);
			
		// Read global values
		double ws						= values.get(0);
		double c						= values.get(1);
		double ksatMult					= values.get(2);
		double bubble					= values.get(3);
		double bulkMult					= values.get(4);
		double soil_bulk				= values.get(5);
		double wrcMult					= values.get(6);
		double rough					= values.get(7);
		double snow_rough				= values.get(8);
		double resid_moist				= values.get(9);
		double cellX					= values.get(10);
		double channelX					= values.get(11);
		double channelKMult				= values.get(12);
		
		resid_moist						= Math.min((1 - 1/soil_bulk)*0.99, resid_moist);
		
		// Read cell variables
		ArrayList<Double> infilt		= new ArrayList<>();
		ArrayList<Double> ds			= new ArrayList<>();
		ArrayList<Double> dsmax			= new ArrayList<>();
		ArrayList<Double> expt			= new ArrayList<>();
		ArrayList<Double> exptMult		= new ArrayList<>();
		ArrayList<Double> ksat			= new ArrayList<>();
		ArrayList<Double> bulk_density	= new ArrayList<>();
		ArrayList<Double> wrc_fract		= new ArrayList<>();
		ArrayList<Double> wpwp_wrc		= new ArrayList<>();
		ArrayList<Double> routingK		= new ArrayList<>();
		ArrayList<Double> directFract	= new ArrayList<>();
		for (int s = 0; s < origSoils.size(); s++)
		{
			infilt.add(			values.get(13 + s*11 + 0));
			ds.add(				values.get(13 + s*11 + 1));
			dsmax.add(			values.get(13 + s*11 + 2));
			expt.add(			values.get(13 + s*11 + 3));
			exptMult.add(		values.get(13 + s*11 + 4));
			ksat.add(			values.get(13 + s*11 + 5));
			bulk_density.add(	values.get(13 + s*11 + 6));
			wrc_fract.add(		values.get(13 + s*11 + 7));
			wpwp_wrc.add(		values.get(13 + s*11 + 8));
			routingK.add(		values.get(13 + s*11 + 9));
			directFract.add(	values.get(13 + s*11 + 10));
		}
		
		// Assign values
		Hashtable<String, MuskingumElement> elements = network.getElements();
		directFractions					= new Hashtable<>();
		for (int s = 0; s < soils.size(); s++)
		{
			Soil soil					= soils.get(s);
			soil.infilt					= infilt.get(s);
			soil.ds						= ds.get(s);
			soil.dsmax					= dsmax.get(s);
			soil.ws						= ws;
			soil.c						= c;
			soil.expt[1]				= expt.get(s);
			soil.expt[0]				= soil.expt[1]*(1 - exptMult.get(s));
			soil.expt[2]				= soil.expt[1]*(1 + exptMult.get(s));
			soil.ksat[1]				= ksat.get(s);
			soil.ksat[0]				= soil.ksat[1]*(1 - ksatMult);
			soil.ksat[2]				= soil.ksat[1]*(1 + ksatMult);
			soil.bubble[0]				= bubble;
			soil.bubble[1]				= bubble;
			soil.bubble[2]				= bubble;
			soil.bulk_density[1]		= bulk_density.get(s);
			soil.bulk_density[0]		= soil.bulk_density[1]*(1 - bulkMult);
			soil.bulk_density[2]		= soil.bulk_density[1]*(1 + bulkMult);
			soil.soil_density[0]		= soil.bulk_density[0]*soil_bulk;
			soil.soil_density[1]		= soil.bulk_density[1]*soil_bulk;
			soil.soil_density[2]		= soil.bulk_density[2]*soil_bulk;
			soil.wrc_fract[1]			= Math.max(resid_moist, wrc_fract.get(s));
			soil.wrc_fract[0]			= Math.max(resid_moist, 
													Math.min(soil.wrc_fract[1]*(1 - wrcMult), 1));
			soil.wrc_fract[2]			= Math.max(resid_moist, 
													Math.min(soil.wrc_fract[1]*(1 + wrcMult), 1));
			soil.wpwp_fract[0]			= Math.max(resid_moist, soil.wrc_fract[0]*wpwp_wrc.get(s));
			soil.wpwp_fract[1]			= Math.max(resid_moist, soil.wrc_fract[1]*wpwp_wrc.get(s));
			soil.wpwp_fract[2]			= Math.max(resid_moist, soil.wrc_fract[2]*wpwp_wrc.get(s));
			soil.rough					= rough;
			soil.snow_rough				= snow_rough;
			soil.resid_moist[0]			= resid_moist;
			soil.resid_moist[1]			= resid_moist;
			soil.resid_moist[2]			= resid_moist;
			
			String routingId			= soil.lat + ", " + soil.lon;
			MuskingumElement element	= elements.get(routingId);
			element.setX(cellX);
			element.setK(routingK.get(s));
			directFractions.put(routingId, directFract.get(s));
		}
		for (String id : elements.keySet())
		{
			MuskingumElement element	= elements.get(id);
			String[] tokens				= id.trim().split(", ");
			if (tokens.length == 1)
			{
				element.setX(channelX);
				element.setK(element.getK()*channelKMult);
			}
		}
		
		return new ModelConfiguration(state, soils, network, directFractions);
	}
	
	public ArrayList<Double> toArray(ModelConfiguration configuration)
	{
		State state					= configuration.state;
		ArrayList<Soil> soils		= configuration.soils;
		MuskingumNetwork network	= configuration.network;
		Hashtable<String, Double> directFractions = configuration.directFractions;
		
		ArrayList<Double> values	= new ArrayList<>();
		
		// Extract global values
		Soil soil					= soils.get(0);
		values.add(soil.ws);
		values.add(soil.c);
		double ksat0				= soil.ksat[0];
		double ksat2				= soil.ksat[2];
		values.add((ksat2 - ksat0)/(ksat0 + ksat2));
		values.add(soil.bubble[1]);
		double bulk0				= soil.bulk_density[0];
		double bulk2				= soil.bulk_density[2];
		values.add((bulk2 - bulk0)/(bulk0 + bulk2));
		values.add(soil.soil_density[1]/soil.bulk_density[1]);
		double wrc0					= soil.wrc_fract[0];
		double wrc2					= soil.wrc_fract[2];
		values.add((wrc2 - wrc0)/(wrc2 + wrc0));
		values.add(soil.rough);
		values.add(soil.snow_rough);
		values.add(soil.resid_moist[1]);
		boolean okCell				= false;
		boolean okChannel			= false;
		Enumeration<String> ids		= network.getElements().keys();
		String channelId			= null;
		while (!(okCell && okChannel))
		{
			String id				= ids.nextElement();
			String[] tokens			= id.trim().split(", ");
			if (!okCell && tokens.length > 1)
			{
				values.add(network.getElements().get(id).getX());
				okCell				= true;
			}
			if (!okChannel && tokens.length == 1)
			{
				channelId			= id;
				okChannel			= true;
			}
		}
		values.add(network.getElements().get(channelId).getX());
		double origK				= origNetwork.getElements().get(channelId).getK();
		double candK				= network.getElements().get(channelId).getK();
		values.add(candK/origK);
		
		// Extract individual cell values
		for (int s = 0; s < soils.size(); s++)
		{
			soil					= soils.get(s);
			values.add(soil.infilt);
			values.add(soil.ds);
			values.add(soil.dsmax);
			values.add(soil.expt[1]);
			double expt0			= soil.expt[0];
			double expt2			= soil.expt[2];
			double exptMult			= expt0 + expt2 == 0.0 ? 0.0 : (expt2 - expt0)/(expt0 + expt2);
			values.add(exptMult);
			values.add(soil.ksat[1]);
			values.add(soil.bulk_density[1]);
			values.add(soil.wrc_fract[1]);
			values.add(soil.wpwp_fract[1]/soil.wrc_fract[1]);
			String id				= soil.lat + ", " + soil.lon;
			values.add(network.getElements().get(id).getK());
			values.add(directFractions.get(id));
		}
		
		// Add state values and return
		values.addAll(state.toArray());
		
		// Verify values
		for (int v = 0; v < values.size(); v++)
			if (Double.isNaN(values.get(v)))
				throw new IllegalArgumentException("Value " + v + " is NaN");
		
		return values;
	}

	@Override
	public NonParametric createDistribution(ArrayList<ContMultiSample> samples)
	{
		NonParametric dist	= null;
		if (distType == OPTIMISTS.TYPE_D_KERNEL || distType == OPTIMISTS.TYPE_F_KERNEL)
		{
			dist			= new MultiVarKernelDensity();
			dist.setWeighted(true);
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
			dist			= new KD_GGMLite(true);
			dist.setSamples(samples);
			if (ggmCreator == null)
				((KD_GGMLite)dist).computeGaussianBW(scaling);
			else
				((KD_GGMLite)dist).computeGaussianBW(scaling, ggmCreator);
		}
		return dist;
	}

}
