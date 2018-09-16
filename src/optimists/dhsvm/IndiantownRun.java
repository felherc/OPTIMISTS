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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.time.Duration;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Scanner;

import dhsvm.MetStation;
import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.Input;
import dhsvm.grid.State;
import dhsvm.stream.Stream;
import dhsvm.stream.StreamCell;
import dhsvm.stream.StreamClass;
import dhsvm.stream.StreamNetwork;
import maestro_mo.MAESTRO;
import maestro_mo.gen.GA;
import maestro_mo.gen.MetroACO;
import optimists.OPTIMISTS;
import probDist.multiVar.tools.GGMLiteCreator;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class IndiantownRun
{

	public static void main(String[] args) throws IOException
	{
		int scenario					= 1;
		int runIndex					= 1;
		String problemName = "Indiantown Run 100-m model data assimilation, scenario " + scenario;
		String outputFolder				= "data/Tests/Scenario " + scenario + "/" + runIndex;
		String modelsFolder				= outputFolder + "/Models";
		String forecastFolder			= outputFolder;
		String inputDataFolder			= "data/Indiantown Run";
		boolean keepModelFiles			= false;
		boolean assimilatorReports		= true;
		boolean maestroReports			= false;
		String hallOfFameFolder			= "";
		boolean computePerformance		= true;
		
		Duration modelTimeStep			= Duration.ofHours(1);
		Duration daTimeStep				= Duration.ofHours(2*7*24);
		//int timeLimit					= 12*60*60*1000;
		int timeLimit					= 15*60*1000;
		
		boolean objQNSE					= false;
		boolean objQMAE					= true;
		boolean objQMARE				= false;
		boolean objIndeppdf				= false;
		boolean objpdf					= false;
		boolean objLogpdf				= false;
		boolean objMDist				= false;
		boolean objMForce				= false;
		boolean objMeanMForce			= true;
		boolean obj2TermCost			= false;
		
		double obsError					= 50;
		double bkgrMultiplier			= 1E-3;
		
		int ensembleSize				= 40;
		int candidateCount				= 200;
		int populationSize				= 25;
		double samplePercentage			= 1.0;
		double rootPercentage			= 0.95;
		double particleGreed			= 0.75;
		int distributionType			= OPTIMISTS.TYPE_GGM_LITE;
		double scaling					= 0.01;
		boolean silverman				= true;
		GGMLiteCreator ggmCreator		= null;
		int threadCount					= 1;
		
		MAESTRO maestro					= new MAESTRO("", 0, null, null, false, true);
		maestro.addGenerator(			new GA());
		maestro.addGenerator(			new MetroACO());
		
		int rows						= 55;
		int cols						= 60;
		int layers						= 3;
		Input input						= defaultInput(inputDataFolder, rows, cols);
		ArrayList<Soil> soils			= defaultSoils();
		ArrayList<Vegetation> vegetations = defaultVegetations();
		String surfaceRoutingFile		= inputDataFolder + "/input/surface_routing.txt";
		StreamNetwork network			= defaultStreamNetwork(surfaceRoutingFile);
		
		ArrayList<MetStation> stations	= new ArrayList<>();
		stations.add(new MetStation(	"NLDAS-2 avg", 2102500.0, 1618500.0, 205.0, 
											inputDataFolder + "/met.txt"));
		
		String optionsFile				= inputDataFolder + "/options.txt";
		String areaFile					= inputDataFolder + "/area.txt";
		String constantsFile			= inputDataFolder + "/constants.txt";
		
		String dhsvmExec				= "data/DHSVM/dhsvm312.exe";
		String streamflowFile			= inputDataFolder + "/Observed flow.txt";
		
		LocalDateTime start				= null;
		LocalDateTime end				= null;
		LocalDateTime forecastEnd		= null;
		if (scenario == 1)
		{
			start						= LocalDateTime.of(2009, 7, 26,  0, 0);
			//start						= LocalDateTime.of(2009, 7, 12,  0, 0);
			//end						= LocalDateTime.of(2009, 7, 26,  6, 0);
			end							= LocalDateTime.of(2009, 8,  9,  0, 0);
			forecastEnd					= LocalDateTime.of(2009, 8, 23,  0, 0);
		}
		else if (scenario == 2)
		{
			start						= LocalDateTime.of(2009, 8, 26,  0, 0);
			//start						= LocalDateTime.of(2009, 8, 12,  0, 0);
			//end						= LocalDateTime.of(2009, 8, 26,  6, 0);
			end							= LocalDateTime.of(2009, 9,  9,  0, 0);
			forecastEnd					= LocalDateTime.of(2009, 9, 23,  0, 0);
		}
		
		ArrayList<State> initialStates	= new ArrayList<>();
		String statesFolder				= inputDataFolder + "/states/scenario " + scenario + "/";
		initialStates.add(new State(start, statesFolder + "state01", rows, cols));
		initialStates.add(new State(start, statesFolder + "state02", rows, cols));
		initialStates.add(new State(start, statesFolder + "state03", rows, cols));
		initialStates.add(new State(start, statesFolder + "state04", rows, cols));
		initialStates.add(new State(start, statesFolder + "state05", rows, cols));
		initialStates.add(new State(start, statesFolder + "state06", rows, cols));
		initialStates.add(new State(start, statesFolder + "state07", rows, cols));
		
		Hashtable<LocalDateTime, Double> streamflowObs = loadDischargeObs(streamflowFile);
		
		DHSVMAssimilator assimilator	= new DHSVMAssimilator(modelTimeStep, input, layers, soils,
				vegetations, network, stations, optionsFile, areaFile, constantsFile, dhsvmExec, 
				objQNSE, objQMAE, objQMARE, objIndeppdf, objpdf, objLogpdf, objMDist, objMForce, 
				objMeanMForce, obj2TermCost, streamflowObs, obsError, bkgrMultiplier);
		assimilator.assimilate(problemName, runIndex, outputFolder, modelsFolder, keepModelFiles, 
				start, maestro, end, daTimeStep, initialStates, ensembleSize, candidateCount, 
				populationSize, samplePercentage, rootPercentage, distributionType, scaling,
				silverman, ggmCreator, particleGreed,  threadCount, assimilatorReports,
				maestroReports, hallOfFameFolder, timeLimit);
		assimilator.forecast(forecastEnd, forecastFolder, threadCount, computePerformance);
	}
	
	private static Hashtable<LocalDateTime, Double> loadDischargeObs(String dischargeFile) 
			throws FileNotFoundException
	{
		Scanner scanner				= new Scanner(new FileInputStream(new File(dischargeFile)));
		Hashtable<LocalDateTime, Double> dischargeObs = new Hashtable<>();
		while (scanner.hasNextLine())
		{
			String line				= scanner.nextLine();
			String[] tokens			= line.split("\t");
			DateTimeFormatter formatter = DateTimeFormatter.ofPattern("MM-dd-yyyy-HH.mm.ss");
			LocalDateTime timeStamp	= LocalDateTime.parse(tokens[0], formatter);
			double value			= Double.valueOf(tokens[1]);
			dischargeObs.put(timeStamp, value);
		}
		scanner.close();
		return dischargeObs;
	}
	
	private static Input defaultInput(String folder, int rows, int cols) 
	{
		Input input	= new Input(rows, cols);
		try 
		{
			input.loadElevationFromFile(	folder + "/input/DEM.bin"			);
			input.loadSoilFromFile(			folder + "/input/Soil.bin"			);
			input.loadSoilDepthFromFile(	folder + "/input/Soil_depth.bin"	);
			input.loadMaskFromFile(			folder + "/input/Mask.bin"			);
			input.loadFlowDirectionFromFile(folder + "/input/Flow_dir.bin"		);
			input.setConstantVegetation(	1									);
		} catch (IOException e) 
		{
			e.printStackTrace();
		}
		return input;
	}
	
	private static ArrayList<Soil> defaultSoils() 
	{
		double[] porosity			= {0.49512301034386846, 0.5914790444942075, 0.6};
		double[] poreSizeDist		= {0.21, 0.21, 0.21};
		double[] bubblingPress		= {0.15, 0.15, 0.15};
		double[] fieldCapacity		= {0.28225696872849537, 0.33718708013472015, 0.39211719154094493};
		double[] wiltingPoint		= {0.11831553195982937, 0.14134095230965138, 0.1643663726594734};
		double[] bulkDensity		= {1569, 1569, 1569};
		double[] vertConductivity	= {0.016616950193884712, 0.016616950193884712, 0.016616950193884712};
		double[] thermalConduct		= {7.114, 6.923, 7.0};
		double[] thermalCapacity	= {1.4E6, 1.4E6, 1.4E6};
		Soil soil1					= new Soil("Sandy loam");
		soil1.setIndex(				1);
		soil1.setLatConductivity(	0.03520397171971311);
		soil1.setExponDecrease(		7.621044496482549);
		soil1.setMaxInfiltration(	3E-5);
		soil1.setDepthThreshold(	0.5);
		soil1.setCapillaryDrive(	0.07);
		soil1.setSurfaceAlbedo(		0.1);
		soil1.setnManning(			0.08);
		soil1.setPorosity(			porosity);
		soil1.setPoreSizeDist(		poreSizeDist);
		soil1.setBubblingPress(		bubblingPress);
		soil1.setFieldCapacity(		fieldCapacity);
		soil1.setWiltingPoint(		wiltingPoint);
		soil1.setBulkDensity(		bulkDensity);
		soil1.setVertConductivity(	vertConductivity);
		soil1.setThermalConduct(	thermalConduct);
		soil1.setThermalCapacity(	thermalCapacity);
		
		double[] porosity2			= {0.21490450619745924, 0.4143730476640377, 0.60};
		double[] poreSizeDist2		= {0.26, 0.26, 0.26};
		double[] bubblingPress2		= {0.21, 0.21, 0.21};
		double[] fieldCapacity2		= {0.14126070755262968, 0.27237506992979726, 0.4034894323069648};
		double[] wiltingPoint2		= {0.04779663916800819, 0.09216018496116592, 0.13652373075432364};
		double[] bulkDensity2		= {1419, 1419, 1419};
		double[] vertConductivity2	= {0.02119972779575767, 0.02119972779575767, 0.02119972779575767};
		double[] thermalConduct2	= {7.114, 6.923, 7.0};
		double[] thermalCapacity2	= {1.4E6, 1.4E6, 1.4E6};
		Soil soil2					= new Soil("Silty loam");
		soil2.setIndex(				2);
		soil2.setLatConductivity(	0.06655098791410533);
		soil2.setExponDecrease(		7.671494624178811);
		soil2.setMaxInfiltration(	3E-5);
		soil2.setDepthThreshold(	0.5);
		soil2.setCapillaryDrive(	0.07);
		soil2.setSurfaceAlbedo(		0.1);
		soil2.setnManning(			0.08);
		soil2.setPorosity(			porosity2);
		soil2.setPoreSizeDist(		poreSizeDist2);
		soil2.setBubblingPress(		bubblingPress2);
		soil2.setFieldCapacity(		fieldCapacity2);
		soil2.setWiltingPoint(		wiltingPoint2);
		soil2.setBulkDensity(		bulkDensity2);
		soil2.setVertConductivity(	vertConductivity2);
		soil2.setThermalConduct(	thermalConduct2);
		soil2.setThermalCapacity(	thermalCapacity2);
		
		double[] porosity3			= {0.43, 0.43, 0.43};
		double[] poreSizeDist3		= {0.19, 0.19, 0.19};
		double[] bubblingPress3		= {0.11, 0.11, 0.11};
		double[] fieldCapacity3		= {0.29, 0.29, 0.29};
		double[] wiltingPoint3		= {0.14, 0.14, 0.14};
		double[] bulkDensity3		= {1485, 1485, 1485};
		double[] vertConductivity3	= {0.01,  0.01,  0.01};
		double[] thermalConduct3	= {7.114, 6.923, 7.0};
		double[] thermalCapacity3	= {1.4E6, 1.4E6, 1.4E6};
		Soil soil3					= new Soil("Loam");
		soil3.setIndex(				3);
		soil3.setLatConductivity(	0.01);
		soil3.setExponDecrease(		3.0);
		soil3.setMaxInfiltration(	1E-5);
		soil3.setDepthThreshold(	0.5);
		soil3.setCapillaryDrive(	0.07);
		soil3.setSurfaceAlbedo(		0.1);
		soil3.setnManning(			0.08);
		soil3.setPorosity(			porosity3);
		soil3.setPoreSizeDist(		poreSizeDist3);
		soil3.setBubblingPress(		bubblingPress3);
		soil3.setFieldCapacity(		fieldCapacity3);
		soil3.setWiltingPoint(		wiltingPoint3);
		soil3.setBulkDensity(		bulkDensity3);
		soil3.setVertConductivity(	vertConductivity3);
		soil3.setThermalConduct(	thermalConduct3);
		soil3.setThermalCapacity(	thermalCapacity3);

		ArrayList<Soil> soils		= new ArrayList<>();
		soils.add(					soil1);
		soils.add(					soil2);
		soils.add(					soil3);
		return soils;
	}

	private static ArrayList<Vegetation> defaultVegetations() 
	{
		ArrayList<Vegetation> vegetations = new ArrayList<>();
		
		Vegetation vegetation		= new Vegetation("Water");
		double[] rootZoneDepth1		= {0.1, 0.25, 0.4};
		double[] usRootFrac1		= {0.0, 0.0, 0.0};
		vegetation.index			= 1;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= false;
		vegetation.impFraction		= 0.0;
		vegetation.rootZoneDepth	= rootZoneDepth1;
		vegetation.usRootFrac		= usRootFrac1;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Urban - Open space");
		double[] height2			= {0.2};
		double[] maxResistance2		= {3000.0};
		double[] minResistance2		= { 120.0};
		double[] moistThres2		= {0.33};
		double[] vaporPressDef2		= {4000.0};
		double[] rpc2				= {0.108};
		double[] rootZoneDepth2		= {0.1, 0.25, 0.4};
		double[] usRootFrac2		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI2		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb2		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 2;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.1;
		vegetation.detenFraction	= 0.0;
		vegetation.detenDecay		= 0.0;
		vegetation.height			= height2;
		vegetation.maxResistance	= maxResistance2;
		vegetation.minResistance	= minResistance2;
		vegetation.moistThres		= moistThres2;
		vegetation.vaporPressDef	= vaporPressDef2;
		vegetation.rpc				= rpc2;
		vegetation.rootZoneDepth	= rootZoneDepth2;
		vegetation.usRootFrac		= usRootFrac2;
		vegetation.usMonthlyLAI		= usMonthlyLAI2;
		vegetation.usMonthlyAlb		= usMonthlyAlb2;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Urban - Low intensity");
		double[] height3			= {0.2};
		double[] maxResistance3		= {3000.0};
		double[] minResistance3		= { 120.0};
		double[] moistThres3		= {   0.33};
		double[] vaporPressDef3		= {4000.0};
		double[] rpc3				= {   0.108};
		double[] rootZoneDepth3		= {0.1, 0.25, 0.4};
		double[] usRootFrac3		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI3		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb3		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 3;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.35;
		vegetation.detenFraction	= 0.0;
		vegetation.detenDecay		= 0.0;
		vegetation.height			= height3;
		vegetation.maxResistance	= maxResistance3;
		vegetation.minResistance	= minResistance3;
		vegetation.moistThres		= moistThres3;
		vegetation.vaporPressDef	= vaporPressDef3;
		vegetation.rpc				= rpc3;
		vegetation.rootZoneDepth	= rootZoneDepth3;
		vegetation.usRootFrac		= usRootFrac3;
		vegetation.usMonthlyLAI		= usMonthlyLAI3;
		vegetation.usMonthlyAlb		= usMonthlyAlb3;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Urban - Medium intensity");
		double[] height4			= {   0.2};
		double[] maxResistance4		= {3000.0};
		double[] minResistance4		= { 120.0};
		double[] moistThres4		= {   0.33};
		double[] vaporPressDef4		= {4000.0};
		double[] rpc4				= {   0.108};
		double[] rootZoneDepth4		= {0.1, 0.25, 0.4};
		double[] usRootFrac4		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI4		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb4		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 4;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.65;
		vegetation.detenFraction	= 0.0;
		vegetation.detenDecay		= 0.0;
		vegetation.height			= height4;
		vegetation.maxResistance	= maxResistance4;
		vegetation.minResistance	= minResistance4;
		vegetation.moistThres		= moistThres4;
		vegetation.vaporPressDef	= vaporPressDef4;
		vegetation.rpc				= rpc4;
		vegetation.rootZoneDepth	= rootZoneDepth4;
		vegetation.usRootFrac		= usRootFrac4;
		vegetation.usMonthlyLAI		= usMonthlyLAI4;
		vegetation.usMonthlyAlb		= usMonthlyAlb4;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Urban - High intensity");
		double[] height5			= {   0.2};
		double[] maxResistance5		= {3000.0};
		double[] minResistance5		= { 120.0};
		double[] moistThres5		= {   0.33};
		double[] vaporPressDef5		= {4000.0};
		double[] rpc5				= {   0.108};
		double[] rootZoneDepth5		= {0.1, 0.25, 0.4};
		double[] usRootFrac5		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI5		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb5		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 5;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.9;
		vegetation.detenFraction	= 0.0;
		vegetation.detenDecay		= 0.0;
		vegetation.height			= height5;
		vegetation.maxResistance	= maxResistance5;
		vegetation.minResistance	= minResistance5;
		vegetation.moistThres		= moistThres5;
		vegetation.vaporPressDef	= vaporPressDef5;
		vegetation.rpc				= rpc5;
		vegetation.rootZoneDepth	= rootZoneDepth5;
		vegetation.usRootFrac		= usRootFrac5;
		vegetation.usMonthlyLAI		= usMonthlyLAI5;
		vegetation.usMonthlyAlb		= usMonthlyAlb5;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Deciduous forest");
		double[] height6			= {30.0, 0.5};
		double[] maxResistance6		= {5000.0, 3000.0};
		double[] minResistance6		= {47.77370440154576, 1514.2638141076832};
		double[] moistThres6		= {0.33, 0.13};
		double[] vaporPressDef6		= {4000.0, 4000.0};
		double[] rpc6				= {0.108, 0.108};
		double[] rootZoneDepth6		= {0.1, 0.25, 0.4};
		double[] osRootFrac6		= {0.2, 0.4, 0.4};
		double[] usRootFrac6		= {0.4, 0.6, 0.0};
		double[] osMonthlyLAI6		= {0.5, 0.5, 0.7, 0.85, 2.5, 3.5, 4.5, 5, 3.5, 1.2, 0.4, 0.3};
		double[] usMonthlyLAI6		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] osMonthlyAlb6		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		double[] usMonthlyAlb6		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 6;
		vegetation.hasOverstory		= true;
		vegetation.hasUnderstory	= true;
		vegetation.fracCover		= 0.6132242485297067;
		vegetation.trunkSpace		= 0.5;
		vegetation.aeroAtten		= 1.5;
		vegetation.radiatAtten		= 0.2;
		vegetation.maxSnowInt		= 0.003;
		vegetation.massRelease		= 0.4;
		vegetation.snowIntEff		= 0.6;
		vegetation.impFraction		= 0.05;
		vegetation.height			= height6;
		vegetation.maxResistance	= maxResistance6;
		vegetation.minResistance	= minResistance6;
		vegetation.moistThres		= moistThres6;
		vegetation.vaporPressDef	= vaporPressDef6;
		vegetation.rpc				= rpc6;
		vegetation.rootZoneDepth	= rootZoneDepth6;
		vegetation.osRootFrac		= osRootFrac6;
		vegetation.usRootFrac		= usRootFrac6;
		vegetation.osMonthlyLAI		= osMonthlyLAI6;
		vegetation.usMonthlyLAI		= usMonthlyLAI6;
		vegetation.osMonthlyAlb		= osMonthlyAlb6;
		vegetation.usMonthlyAlb		= usMonthlyAlb6;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Evergreen forest");
		double[] height7			= {30.0, 0.5};
		double[] maxResistance7		= {5000.0, 3000.0};
		double[] minResistance7		= {666.6, 666.6};
		double[] moistThres7		= {0.33, 0.13};
		double[] vaporPressDef7		= {4000.0, 4000.0};
		double[] rpc7				= {0.108, 0.108};
		double[] rootZoneDepth7		= {0.1, 0.25, 0.4};
		double[] osRootFrac7		= {0.2, 0.4, 0.4};
		double[] usRootFrac7		= {0.4, 0.6, 0.0};
		double[] osMonthlyLAI7		= {4, 4, 4, 4, 4, 4, 4.5, 5, 4.5, 4, 4, 4};
		double[] usMonthlyLAI7		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] osMonthlyAlb7		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		double[] usMonthlyAlb7		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 7;
		vegetation.hasOverstory		= true;
		vegetation.hasUnderstory	= true;
		vegetation.fracCover		= 0.9;
		vegetation.trunkSpace		= 0.5;
		vegetation.aeroAtten		= 1.5;
		vegetation.radiatAtten		= 0.2;
		vegetation.maxSnowInt		= 0.003;
		vegetation.massRelease		= 0.4;
		vegetation.snowIntEff		= 0.6;
		vegetation.impFraction		= 0.05;
		vegetation.height			= height7;
		vegetation.maxResistance	= maxResistance7;
		vegetation.minResistance	= minResistance7;
		vegetation.moistThres		= moistThres7;
		vegetation.vaporPressDef	= vaporPressDef7;
		vegetation.rpc				= rpc7;
		vegetation.rootZoneDepth	= rootZoneDepth7;
		vegetation.osRootFrac		= osRootFrac7;
		vegetation.usRootFrac		= usRootFrac7;
		vegetation.osMonthlyLAI		= osMonthlyLAI7;
		vegetation.usMonthlyLAI		= usMonthlyLAI7;
		vegetation.osMonthlyAlb		= osMonthlyAlb7;
		vegetation.usMonthlyAlb		= usMonthlyAlb7;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Mixed forest");
		double[] height8			= {20.0, 0.5};
		double[] maxResistance8		= {5000.0, 600.0};
		double[] minResistance8		= {200.0, 200.0};
		double[] moistThres8		= {0.33, 0.13};
		double[] vaporPressDef8		= {4000.0, 4000.0};
		double[] rpc8				= {0.108, 0.108};
		double[] rootZoneDepth8		= {0.1, 0.25, 0.4};
		double[] osRootFrac8		= {0.2, 0.4, 0.4};
		double[] usRootFrac8		= {0.4, 0.6, 0.0};
		double[] osMonthlyLAI8		= {0.5, 0.5, 0.7, 0.85, 2.5, 3.5, 4.5, 5, 3.5, 1.2, 0.4, 0.3};
		double[] usMonthlyLAI8		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] osMonthlyAlb8		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		double[] usMonthlyAlb8		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 8;
		vegetation.hasOverstory		= true;
		vegetation.hasUnderstory	= true;
		vegetation.fracCover		= 0.8;
		vegetation.trunkSpace		= 0.4;
		vegetation.aeroAtten		= 0.5;
		vegetation.radiatAtten		= 0.2;
		vegetation.maxSnowInt		= 0.003;
		vegetation.massRelease		= 0.4;
		vegetation.snowIntEff		= 0.6;
		vegetation.impFraction		= 0.05;
		vegetation.height			= height8;
		vegetation.maxResistance	= maxResistance8;
		vegetation.minResistance	= minResistance8;
		vegetation.moistThres		= moistThres8;
		vegetation.vaporPressDef	= vaporPressDef8;
		vegetation.rpc				= rpc8;
		vegetation.rootZoneDepth	= rootZoneDepth8;
		vegetation.osRootFrac		= osRootFrac8;
		vegetation.usRootFrac		= usRootFrac8;
		vegetation.osMonthlyLAI		= osMonthlyLAI8;
		vegetation.usMonthlyLAI		= usMonthlyLAI8;
		vegetation.osMonthlyAlb		= osMonthlyAlb8;
		vegetation.usMonthlyAlb		= usMonthlyAlb8;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Pasture/Hay");
		double[] height9			= {   0.5};
		double[] maxResistance9		= { 600.0};
		double[] minResistance9		= { 200.0};
		double[] moistThres9		= {   0.33};
		double[] vaporPressDef9		= {4000.0};
		double[] rpc9				= {   0.108};
		double[] rootZoneDepth9		= {0.1, 0.25, 0.4};
		double[] usRootFrac9		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI9		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb9		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 9;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.05;
		vegetation.height			= height9;
		vegetation.maxResistance	= maxResistance9;
		vegetation.minResistance	= minResistance9;
		vegetation.moistThres		= moistThres9;
		vegetation.vaporPressDef	= vaporPressDef9;
		vegetation.rpc				= rpc9;
		vegetation.rootZoneDepth	= rootZoneDepth9;
		vegetation.usRootFrac		= usRootFrac9;
		vegetation.usMonthlyLAI		= usMonthlyLAI9;
		vegetation.usMonthlyAlb		= usMonthlyAlb9;
		vegetations.add(			vegetation);
		
		vegetation					= new Vegetation("Cultivated crops");
		double[] height10			= {   1.0};
		double[] maxResistance10	= { 600.0};
		double[] minResistance10	= { 120.0};
		double[] moistThres10		= {   0.33};
		double[] vaporPressDef10	= {4000.0};
		double[] rpc10				= {   0.108};
		double[] rootZoneDepth10	= {0.1, 0.25, 0.4};
		double[] usRootFrac10		= {0.4, 0.6, 0.0};
		double[] usMonthlyLAI10		= {0.162, 0.150, 0.162, 0.275, 0.525, 1.037, 0.962, 0.550, 0.338, 0.200, 0.138, 0.138};
		double[] usMonthlyAlb10		= {0.05, 0.05, 0.05, 0.07, 0.04, 0.04, 0.04, 0.04, 0.04, 0.05, 0.05, 0.05};
		vegetation.index			= 10;
		vegetation.hasOverstory		= false;
		vegetation.hasUnderstory	= true;
		vegetation.impFraction		= 0.05;
		vegetation.height			= height10;
		vegetation.maxResistance	= maxResistance10;
		vegetation.minResistance	= minResistance10;
		vegetation.moistThres		= moistThres10;
		vegetation.vaporPressDef	= vaporPressDef10;
		vegetation.rpc				= rpc10;
		vegetation.rootZoneDepth	= rootZoneDepth10;
		vegetation.usRootFrac		= usRootFrac10;
		vegetation.usMonthlyLAI		= usMonthlyLAI10;
		vegetation.usMonthlyAlb		= usMonthlyAlb10;
		vegetations.add(			vegetation);
		
		return vegetations;
	}

	private static StreamNetwork defaultStreamNetwork(String surfaceRoutingFile)
	{
		StreamNetwork network		= new StreamNetwork();
		StreamClass class1   	= new StreamClass(1, 1.84734912505056, 0.277102368757584, 0.157033709931452, 0.00001);
		StreamClass class2   	= new StreamClass(2, 2.02507410365325, 0.303761115547988, 0.157033709931452, 0.00001);
		StreamClass class3   	= new StreamClass(3, 1.94605296805197, 0.291907945207796, 0.157033709931452, 0.00001);
		StreamClass class4   	= new StreamClass(4, 8.90609792893338, 1.33591468934001, 0.157033709931452, 0.00001);
		StreamClass class5   	= new StreamClass(5, 8.18225247031899, 1.22733787054785, 0.157033709931452, 0.00001);
		StreamClass class6  	= new StreamClass(6, 8.55359094484427, 1.28303864172664, 0.157033709931452, 0.00001);
		StreamClass class7  	= new StreamClass(7, 8.59566088426779, 1.28934913264017, 0.157033709931452, 0.00001);
		StreamClass class8   	= new StreamClass(8, 5.61934156378601, 0.842901234567901, 0.157033709931452, 0.00001);
		StreamClass class9   	= new StreamClass(9, 7.31292478959786, 1.09693871843968, 0.157033709931452, 0.00001);
		StreamClass class10   	= new StreamClass(10, 4.60428176231937, 0.690642264347906, 0.157033709931452, 0.00001);
		StreamClass class11   	= new StreamClass(11, 3.62111231932409, 0.543166847898613, 0.157033709931452, 0.00001);
		StreamClass class12   	= new StreamClass(12, 1.71886173295642, 0.257829259943463, 0.157033709931452, 0.00001);
		StreamClass class13   	= new StreamClass(13, 1.95161083856062, 0.292741625784094, 0.157033709931452, 0.00001);
		StreamClass class14   	= new StreamClass(14, 3.2201070562821, 0.483016058442316, 0.157033709931452, 0.00001);
		StreamClass class15   	= new StreamClass(15, 2.13446317052877, 0.320169475579316, 0.157033709931452, 0.00001);
		StreamClass class16   	= new StreamClass(16, 2.10150574184511, 0.315225861276766, 0.157033709931452, 0.00001);
		StreamClass class17   	= new StreamClass(17, 1.66671792444021, 0.250007688666032, 0.157033709931452, 0.00001);
		StreamClass class18   	= new StreamClass(18, 1.7263781537686, 0.258956723065291, 0.157033709931452, 0.00001);
		StreamClass class19   	= new StreamClass(19, 3.58366435960201, 0.537549653940301, 0.157033709931452, 0.00001);
		StreamClass class20   	= new StreamClass(20, 2.30462045322996, 0.345693067984494, 0.157033709931452, 0.00001);
		StreamClass class21   	= new StreamClass(21, 1.65061978854322, 0.247592968281483, 0.157033709931452, 0.00001);
		Stream stream4   		= new Stream(4, class4, 14, 0.030712, 1086.801355, 0, null);
		Stream stream7   		= new Stream(7, class7, 13, 0.0286, 651.275898, 0, stream4);
		Stream stream5   		= new Stream(5, class5, 12, 0.032344, 550.168818, 0, stream7);
		Stream stream1   		= new Stream(1, class1, 1, 0.06308, 1356.957789, 0, stream5);
		Stream stream6  		= new Stream(6, class6, 11, 0.025933, 77.843528, 0, stream5);
		Stream stream9   		= new Stream(9, class9, 10, 0.025074, 586.38957, 0, stream6);
		Stream stream10   		= new Stream(10, class10, 6, 0.029015, 181.676797, 0, stream9);
		Stream stream2   		= new Stream(2, class2, 1, 0.054738, 1766.883859, 0, stream10);
		Stream stream8  		= new Stream(8, class8, 5, 0.013922, 142.689916, 0, stream10);
		Stream stream11   		= new Stream(11, class11, 4, 0.024637, 482.611251, 0, stream8);
		Stream stream3   		= new Stream(3, class3, 1, 0.026781, 1188.377361, 0, stream11);
		Stream stream14   		= new Stream(14, class14, 3, 0.021299, 850.956593, 0, stream11);
		Stream stream12   		= new Stream(12, class12, 1, 0.051881, 703.085261, 0, stream14);
		Stream stream13   		= new Stream(13, class13, 1, 0.11229, 1154.302779, 0, stream4);
		Stream stream15   		= new Stream(15, class15, 2, 0.039695, 199.727028, 0, stream14);
		Stream stream16   		= new Stream(16, class16, 1, 0.025547, 1771.866, 0, stream6);
		Stream stream17   		= new Stream(17, class17, 1, 0.109595, 697.966617, 0, stream15);
		Stream stream19   		= new Stream(19, class19, 3, 0.016868, 982.879469, 0, stream9);
		Stream stream20   		= new Stream(20, class20, 2, 0.020686, 394.168307, 0, stream19);
		Stream stream21   		= new Stream(21, class21, 1, 0.05044, 368.238951, 0, stream20);
		Stream stream18   		= new Stream(18, class18, 1, 0.029937, 300.998151, 0, stream19);
		network.addStreamClass(	class1);
		network.addStreamClass(	class2);
		network.addStreamClass(	class3);
		network.addStreamClass(	class4);
		network.addStreamClass(	class5);
		network.addStreamClass(	class6);
		network.addStreamClass(	class7);
		network.addStreamClass(	class8);
		network.addStreamClass(	class9);
		network.addStreamClass(	class10);
		network.addStreamClass(	class11);
		network.addStreamClass(	class12);
		network.addStreamClass(	class13);
		network.addStreamClass(	class14);
		network.addStreamClass(	class15);
		network.addStreamClass(	class16);
		network.addStreamClass(	class17);
		network.addStreamClass(	class18);
		network.addStreamClass(	class19);
		network.addStreamClass(	class20);
		network.addStreamClass(	class21);
		network.addStream(		stream1);
		network.addStream(		stream2);
		network.addStream(		stream3);
		network.addStream(		stream4);
		network.addStream(		stream5);
		network.addStream(		stream6);
		network.addStream(		stream7);
		network.addStream(		stream8);
		network.addStream(		stream9);
		network.addStream(		stream10);
		network.addStream(		stream11);
		network.addStream(		stream12);
		network.addStream(		stream13);
		network.addStream(		stream14);
		network.addStream(		stream15);
		network.addStream(		stream16);
		network.addStream(		stream17);
		network.addStream(		stream18);
		network.addStream(		stream19);
		network.addStream(		stream20);
		network.addStream(		stream21);
		network.addStreamCell(  new StreamCell(9, 23, stream17, 24.185549));
		network.addStreamCell(  new StreamCell(9, 24, stream17, 29.648662));
		network.addStreamCell(  new StreamCell(10, 23, stream17, 131.414616));
		network.addStreamCell(  new StreamCell(11, 22, stream17, 53.879806));
		network.addStreamCell(  new StreamCell(11, 23, stream17, 62.883678));
		network.addStreamCell(  new StreamCell(12, 22, stream17, 109.437885));
		network.addStreamCell(  new StreamCell(12, 28, stream12, 54.021471));
		network.addStreamCell(  new StreamCell(12, 29, stream12, 37.1495459999999));
		network.addStreamCell(  new StreamCell(13, 22, stream17, 113.961241));
		network.addStreamCell(  new StreamCell(13, 27, stream12, 97.742612));
		network.addStreamCell(  new StreamCell(13, 28, stream12, 68.607872));
		network.addStreamCell(  new StreamCell(13, 29, stream12, 88.9983030000001));
		network.addStreamCell(  new StreamCell(14, 22, stream17, 64.90881));
		network.addStreamCell(  new StreamCell(14, 23, stream17, 10.003517));
		network.addStreamCell(  new StreamCell(14, 24, stream17, 126.140752));
		network.addStreamCell(  new StreamCell(14, 26, stream12, 117.270868));
		network.addStreamCell(  new StreamCell(14, 27, stream12, 50.77664));
		network.addStreamCell(  new StreamCell(14, 28, stream12, 51.100156));
		network.addStreamCell(  new StreamCell(15, 23, stream17, 115.043262));
		network.addStreamCell(  new StreamCell(15, 23, stream15, 75.129878));
		network.addStreamCell(  new StreamCell(15, 25, stream12, 106.157845));
		network.addStreamCell(  new StreamCell(15, 26, stream12, 106.157845));
		network.addStreamCell(  new StreamCell(15, 27, stream12, 13.202008));
		network.addStreamCell(  new StreamCell(15, 33, stream3, 49.007161));
		network.addStreamCell(  new StreamCell(15, 34, stream3, 75.860204));
		network.addStreamCell(  new StreamCell(16, 23, stream15, 24.367332));
		network.addStreamCell(  new StreamCell(16, 24, stream15, 142.712825));
		network.addStreamCell(  new StreamCell(16, 24, stream12, 115.974734));
		network.addStreamCell(  new StreamCell(16, 24, stream14, 15.558457));
		network.addStreamCell(  new StreamCell(16, 25, stream12, 115.974734));
		network.addStreamCell(  new StreamCell(16, 33, stream3, 137.354588));
		network.addStreamCell(  new StreamCell(17, 24, stream14, 108.192676));
		network.addStreamCell(  new StreamCell(17, 32, stream3, 136.295238));
		network.addStreamCell(  new StreamCell(17, 33, stream3, 3.83601800000002));
		network.addStreamCell(  new StreamCell(17, 34, stream3, 80.8170620000001));
		network.addStreamCell(  new StreamCell(18, 24, stream14, 109.291509));
		network.addStreamCell(  new StreamCell(18, 31, stream3, 103.283718));
		network.addStreamCell(  new StreamCell(18, 32, stream3, 42.787286));
		network.addStreamCell(  new StreamCell(18, 33, stream3, 118.635652));
		network.addStreamCell(  new StreamCell(18, 34, stream3, 45.3556689999998));
		network.addStreamCell(  new StreamCell(18, 38, stream2, 74.975382));
		network.addStreamCell(  new StreamCell(18, 39, stream2, 39.261167));
		network.addStreamCell(  new StreamCell(19, 24, stream14, 60.479933));
		network.addStreamCell(  new StreamCell(19, 25, stream14, 98.834953));
		network.addStreamCell(  new StreamCell(19, 30, stream3, 57.311387));
		network.addStreamCell(  new StreamCell(19, 31, stream3, 88.761103));
		network.addStreamCell(  new StreamCell(19, 32, stream3, 126.172731));
		network.addStreamCell(  new StreamCell(19, 33, stream3, 7.53707899999995));
		network.addStreamCell(  new StreamCell(19, 37, stream2, 89.753485));
		network.addStreamCell(  new StreamCell(19, 38, stream2, 39.261167));
		network.addStreamCell(  new StreamCell(19, 39, stream2, 126.188389));
		network.addStreamCell(  new StreamCell(20, 25, stream14, 39.883564));
		network.addStreamCell(  new StreamCell(20, 26, stream14, 146.924002));
		network.addStreamCell(  new StreamCell(20, 29, stream3, 63.907922));
		network.addStreamCell(  new StreamCell(20, 30, stream3, 89.492525));
		network.addStreamCell(  new StreamCell(20, 31, stream3, 95.89122));
		network.addStreamCell(  new StreamCell(20, 36, stream2, 109.950188));
		network.addStreamCell(  new StreamCell(20, 37, stream2, 56.645699));
		network.addStreamCell(  new StreamCell(20, 38, stream2, 97.5915300000001));
		network.addStreamCell(  new StreamCell(21, 26, stream14, 4.14109300000007));
		network.addStreamCell(  new StreamCell(21, 27, stream14, 127.265234));
		network.addStreamCell(  new StreamCell(21, 29, stream3, 116.767466));
		network.addStreamCell(  new StreamCell(21, 30, stream3, 58.072631));
		network.addStreamCell(  new StreamCell(21, 35, stream2, 133.099251));
		network.addStreamCell(  new StreamCell(21, 36, stream2, 30.813715));
		network.addStreamCell(  new StreamCell(21, 37, stream2, 59.8229629999996));
		network.addStreamCell(  new StreamCell(22, 27, stream3, 102.739901));
		network.addStreamCell(  new StreamCell(22, 27, stream14, 58.174011));
		network.addStreamCell(  new StreamCell(22, 28, stream3, 102.739901));
		network.addStreamCell(  new StreamCell(22, 29, stream3, 50.747963));
		network.addStreamCell(  new StreamCell(22, 34, stream2, 131.914478));
		network.addStreamCell(  new StreamCell(22, 35, stream2, 19.75273));
		network.addStreamCell(  new StreamCell(22, 36, stream2, 22.0543950000001));
		network.addStreamCell(  new StreamCell(23, 27, stream11, 97.438144));
		network.addStreamCell(  new StreamCell(23, 27, stream3, 72.4984449999999));
		network.addStreamCell(  new StreamCell(23, 27, stream14, 58.174011));
		network.addStreamCell(  new StreamCell(23, 28, stream11, 6.33377400000001));
		network.addStreamCell(  new StreamCell(23, 33, stream2, 80.262554));
		network.addStreamCell(  new StreamCell(23, 34, stream2, 51.158968));
		network.addStreamCell(  new StreamCell(24, 28, stream11, 132.221872));
		network.addStreamCell(  new StreamCell(24, 32, stream2, 90.718468));
		network.addStreamCell(  new StreamCell(24, 33, stream2, 90.718468));
		network.addStreamCell(  new StreamCell(25, 28, stream11, 109.30279));
		network.addStreamCell(  new StreamCell(25, 29, stream11, 8.00533999999999));
		network.addStreamCell(  new StreamCell(25, 30, stream2, 2.83148899999969));
		network.addStreamCell(  new StreamCell(25, 31, stream2, 126.188389));
		network.addStreamCell(  new StreamCell(25, 32, stream2, 131.423396));
		network.addStreamCell(  new StreamCell(26, 20, stream18, 10.664469));
		network.addStreamCell(  new StreamCell(26, 28, stream11, 89.028837));
		network.addStreamCell(  new StreamCell(26, 29, stream11, 13.086223));
		network.addStreamCell(  new StreamCell(26, 29, stream2, 40.6000569999999));
		network.addStreamCell(  new StreamCell(26, 30, stream2, 123.3569));
		network.addStreamCell(  new StreamCell(26, 31, stream2, 110.961721));
		network.addStreamCell(  new StreamCell(26, 32, stream2, 20.462714));
		network.addStreamCell(  new StreamCell(27, 11, stream21, 8.061357));
		network.addStreamCell(  new StreamCell(27, 12, stream21, 110.529572));
		network.addStreamCell(  new StreamCell(27, 13, stream21, 102.781013));
		network.addStreamCell(  new StreamCell(27, 14, stream21, 108.618165));
		network.addStreamCell(  new StreamCell(27, 15, stream21, 14.000453));
		network.addStreamCell(  new StreamCell(27, 19, stream18, 139.093295));
		network.addStreamCell(  new StreamCell(27, 20, stream18, 6.991835));
		network.addStreamCell(  new StreamCell(27, 28, stream8, 78.368625));
		network.addStreamCell(  new StreamCell(27, 28, stream2, 74.3988319999999));
		network.addStreamCell(  new StreamCell(27, 28, stream11, 56.889294));
		network.addStreamCell(  new StreamCell(27, 29, stream2, 71.844083));
		network.addStreamCell(  new StreamCell(27, 29, stream11, 13.969037));
		network.addStreamCell(  new StreamCell(27, 30, stream2, 43.8992070000002));
		network.addStreamCell(  new StreamCell(27, 31, stream2, 111.553634));
		network.addStreamCell(  new StreamCell(27, 38, stream1, 82.690265));
		network.addStreamCell(  new StreamCell(28, 15, stream20, 102.711591));
		network.addStreamCell(  new StreamCell(28, 15, stream21, 45.101911));
		network.addStreamCell(  new StreamCell(28, 16, stream20, 102.116337));
		network.addStreamCell(  new StreamCell(28, 17, stream20, 15.864449));
		network.addStreamCell(  new StreamCell(28, 18, stream18, 62.318191));
		network.addStreamCell(  new StreamCell(28, 19, stream18, 38.686659));
		network.addStreamCell(  new StreamCell(28, 20, stream18, 68.9008220000001));
		network.addStreamCell(  new StreamCell(28, 27, stream10, 8.63085500000001));
		network.addStreamCell(  new StreamCell(28, 28, stream10, 75.758794));
		network.addStreamCell(  new StreamCell(28, 28, stream2, 74.3988319999999));
		network.addStreamCell(  new StreamCell(28, 28, stream8, 49.225168));
		network.addStreamCell(  new StreamCell(28, 29, stream2, 38.686659));
		network.addStreamCell(  new StreamCell(28, 30, stream2, 86.1773739999999));
		network.addStreamCell(  new StreamCell(28, 37, stream1, 103.557578));
		network.addStreamCell(  new StreamCell(28, 38, stream1, 103.557578));
		network.addStreamCell(  new StreamCell(29, 17, stream20, 87.183958));
		network.addStreamCell(  new StreamCell(29, 18, stream20, 111.109768));
		network.addStreamCell(  new StreamCell(29, 18, stream18, 59.621897));
		network.addStreamCell(  new StreamCell(29, 18, stream19, 24.750877));
		network.addStreamCell(  new StreamCell(29, 19, stream19, 106.853112));
		network.addStreamCell(  new StreamCell(29, 19, stream18, 31.076415));
		network.addStreamCell(  new StreamCell(29, 20, stream19, 104.09001));
		network.addStreamCell(  new StreamCell(29, 21, stream19, 118.744854));
		network.addStreamCell(  new StreamCell(29, 22, stream19, 105.734451));
		network.addStreamCell(  new StreamCell(29, 23, stream19, 116.44201));
		network.addStreamCell(  new StreamCell(29, 24, stream19, 105.077304));
		network.addStreamCell(  new StreamCell(29, 25, stream19, 107.310765));
		network.addStreamCell(  new StreamCell(29, 26, stream19, 103.103166));
		network.addStreamCell(  new StreamCell(29, 27, stream19, 115.190129));
		network.addStreamCell(  new StreamCell(29, 27, stream10, 108.886789));
		network.addStreamCell(  new StreamCell(29, 27, stream9, 13.894343));
		network.addStreamCell(  new StreamCell(29, 28, stream10, 126.192768));
		network.addStreamCell(  new StreamCell(29, 28, stream9, 3.155288));
		network.addStreamCell(  new StreamCell(29, 36, stream1, 57.369354));
		network.addStreamCell(  new StreamCell(29, 37, stream1, 79.511048));
		network.addStreamCell(  new StreamCell(30, 28, stream9, 128.588511));
		network.addStreamCell(  new StreamCell(30, 35, stream1, 44.661104));
		network.addStreamCell(  new StreamCell(30, 36, stream1, 109.161786));
		network.addStreamCell(  new StreamCell(31, 28, stream9, 60.097152));
		network.addStreamCell(  new StreamCell(31, 29, stream9, 63.726354));
		network.addStreamCell(  new StreamCell(31, 33, stream1, 17.5628220000001));
		network.addStreamCell(  new StreamCell(31, 34, stream1, 67.473938));
		network.addStreamCell(  new StreamCell(31, 35, stream1, 106.581404));
		network.addStreamCell(  new StreamCell(32, 29, stream9, 116.105885));
		network.addStreamCell(  new StreamCell(32, 32, stream1, 55.2552620000001));
		network.addStreamCell(  new StreamCell(32, 33, stream1, 97.588379));
		network.addStreamCell(  new StreamCell(32, 34, stream1, 73.277236));
		network.addStreamCell(  new StreamCell(33, 29, stream9, 113.637049));
		network.addStreamCell(  new StreamCell(33, 31, stream1, 92.9477019999999));
		network.addStreamCell(  new StreamCell(33, 32, stream1, 117.945766));
		network.addStreamCell(  new StreamCell(33, 33, stream1, 49.204968));
		network.addStreamCell(  new StreamCell(34, 25, stream16, 82.922988));
		network.addStreamCell(  new StreamCell(34, 26, stream16, 117.826754));
		network.addStreamCell(  new StreamCell(34, 27, stream16, 105.401687));
		network.addStreamCell(  new StreamCell(34, 28, stream16, 107.583703));
		network.addStreamCell(  new StreamCell(34, 29, stream9, 87.184987));
		network.addStreamCell(  new StreamCell(34, 29, stream16, 129.136163));
		network.addStreamCell(  new StreamCell(34, 29, stream1, 4.42790300000001));
		network.addStreamCell(  new StreamCell(34, 29, stream6, 4.42790300000001));
		network.addStreamCell(  new StreamCell(34, 30, stream1, 42.845811));
		network.addStreamCell(  new StreamCell(34, 31, stream1, 127.804642));
		network.addStreamCell(  new StreamCell(34, 32, stream1, 28.8457079999998));
		network.addStreamCell(  new StreamCell(35, 24, stream16, 56.8318979999999));
		network.addStreamCell(  new StreamCell(35, 25, stream16, 86.915653));
		network.addStreamCell(  new StreamCell(35, 29, stream1, 136.097636));
		network.addStreamCell(  new StreamCell(35, 29, stream6, 120.482133));
		network.addStreamCell(  new StreamCell(35, 29, stream5, 66.472295));
		network.addStreamCell(  new StreamCell(35, 30, stream1, 65.9732949999998));
		network.addStreamCell(  new StreamCell(36, 23, stream16, 127.708125));
		network.addStreamCell(  new StreamCell(36, 24, stream16, 71.82104));
		network.addStreamCell(  new StreamCell(36, 29, stream5, 105.994187));
		network.addStreamCell(  new StreamCell(37, 21, stream16, 4.90806499999997));
		network.addStreamCell(  new StreamCell(37, 22, stream16, 152.744709));
		network.addStreamCell(  new StreamCell(37, 23, stream16, 0.0897650000000567));
		network.addStreamCell(  new StreamCell(37, 28, stream5, 18.090825));
		network.addStreamCell(  new StreamCell(37, 29, stream5, 118.531918));
		network.addStreamCell(  new StreamCell(38, 21, stream16, 25.7333520000001));
		network.addStreamCell(  new StreamCell(38, 22, stream16, 10.013011));
		network.addStreamCell(  new StreamCell(38, 28, stream5, 124.105389));
		network.addStreamCell(  new StreamCell(38, 29, stream5, 73.7120170000001));
		network.addStreamCell(  new StreamCell(39, 19, stream16, 29.718999));
		network.addStreamCell(  new StreamCell(39, 20, stream16, 143.827327));
		network.addStreamCell(  new StreamCell(39, 21, stream16, 10.057856));
		network.addStreamCell(  new StreamCell(39, 28, stream5, 19.605369));
		network.addStreamCell(  new StreamCell(39, 29, stream5, 52.5106049999999));
		network.addStreamCell(  new StreamCell(40, 19, stream16, 29.718999));
		network.addStreamCell(  new StreamCell(40, 20, stream16, 14.285785));
		network.addStreamCell(  new StreamCell(40, 27, stream7, 121.482064));
		network.addStreamCell(  new StreamCell(40, 28, stream5, 19.605369));
		network.addStreamCell(  new StreamCell(40, 28, stream7, 2.42765));
		network.addStreamCell(  new StreamCell(41, 18, stream16, 95.331933));
		network.addStreamCell(  new StreamCell(41, 19, stream16, 70.338507));
		network.addStreamCell(  new StreamCell(41, 26, stream7, 11.635389));
		network.addStreamCell(  new StreamCell(41, 27, stream7, 29.812559));
		network.addStreamCell(  new StreamCell(42, 17, stream16, 77.013993));
		network.addStreamCell(  new StreamCell(42, 18, stream16, 69.364814));
		network.addStreamCell(  new StreamCell(42, 26, stream7, 81.092535));
		network.addStreamCell(  new StreamCell(42, 27, stream7, 0.159464000000014));
		network.addStreamCell(  new StreamCell(43, 17, stream16, 46.113696));
		network.addStreamCell(  new StreamCell(43, 24, stream13, 85.298994));
		network.addStreamCell(  new StreamCell(43, 25, stream13, 107.624416));
		network.addStreamCell(  new StreamCell(43, 26, stream7, 61.097458));
		network.addStreamCell(  new StreamCell(43, 26, stream13, 10.689965));
		network.addStreamCell(  new StreamCell(43, 27, stream7, 59.885059));
		network.addStreamCell(  new StreamCell(44, 21, stream13, 6.46839399999999));
		network.addStreamCell(  new StreamCell(44, 22, stream13, 122.602522));
		network.addStreamCell(  new StreamCell(44, 23, stream13, 107.959979));
		network.addStreamCell(  new StreamCell(44, 24, stream13, 23.514504));
		network.addStreamCell(  new StreamCell(44, 26, stream13, 121.875446));
		network.addStreamCell(  new StreamCell(44, 26, stream7, 20.286176));
		network.addStreamCell(  new StreamCell(44, 27, stream7, 121.14344));
		network.addStreamCell(  new StreamCell(44, 27, stream13, 52.109506));
		network.addStreamCell(  new StreamCell(44, 28, stream7, 85.364978));
		network.addStreamCell(  new StreamCell(44, 28, stream13, 85.3649770000002));
		network.addStreamCell(  new StreamCell(45, 20, stream13, 29.402063));
		network.addStreamCell(  new StreamCell(45, 21, stream13, 140.300346));
		network.addStreamCell(  new StreamCell(45, 27, stream13, 81.3430860000001));
		network.addStreamCell(  new StreamCell(45, 27, stream7, 30.2523150000001));
		network.addStreamCell(  new StreamCell(45, 28, stream4, 122.938455));
		network.addStreamCell(  new StreamCell(45, 28, stream7, 41.16814));
		network.addStreamCell(  new StreamCell(45, 28, stream13, 40.6864859999998));
		network.addStreamCell(  new StreamCell(45, 29, stream4, 43.844187));
		network.addStreamCell(  new StreamCell(46, 19, stream13, 56.772895));
		network.addStreamCell(  new StreamCell(46, 20, stream13, 111.323142));
		network.addStreamCell(  new StreamCell(46, 29, stream4, 102.332712));
		network.addStreamCell(  new StreamCell(46, 30, stream4, 64.322216));
		network.addStreamCell(  new StreamCell(47, 18, stream13, 11.75794));
		network.addStreamCell(  new StreamCell(47, 19, stream13, 83.950734));
		network.addStreamCell(  new StreamCell(47, 30, stream4, 74.112606));
		network.addStreamCell(  new StreamCell(47, 31, stream4, 72.653759));
		network.addStreamCell(  new StreamCell(48, 31, stream4, 79.311391));
		network.addStreamCell(  new StreamCell(48, 32, stream4, 115.95048));
		network.addStreamCell(  new StreamCell(48, 33, stream4, 14.4860609999999));
		network.addStreamCell(  new StreamCell(49, 33, stream4, 116.893818));
		network.addStreamCell(  new StreamCell(49, 34, stream4, 101.046001));
		network.addStreamCell(  new StreamCell(50, 34, stream4, 30.516395));
		network.addStreamCell(  new StreamCell(50, 35, stream4, 119.67886));
		network.addStreamCell(  new StreamCell(51, 35, stream4, 56.897588));
		
		try 
		{
			if (!surfaceRoutingFile.equals(""))
				network.loadSurfaceRouting(surfaceRoutingFile);
		} catch (FileNotFoundException e)
		{
			e.printStackTrace();
		}
		
		return network;
	}

}
