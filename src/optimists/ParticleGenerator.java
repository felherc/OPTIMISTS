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

package optimists;

import java.util.ArrayList;
import java.util.Queue;

import maestro_mo.ContVar;
import maestro_mo.DiscVar;
import maestro_mo.Objective;
import maestro_mo.gen.Generator;
import maestro_mo.pop.Population;
import maestro_mo.solution.SolutionRoot;
import probDist.multiVar.NonParametric;
import utilities.Utilities;

/**
 * Generates solution roots by taking root samples from the current/source state probability
 * distribution and from randomly sampling it
 * 
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class ParticleGenerator implements Generator
{

	// --------------------------------------------------------------------------------------------
	// Constants
	// --------------------------------------------------------------------------------------------
	
	/**
	 * The identifier of the generator
	 */
	public final static String ID = "OPTIMISTS' source distribution sampler";
	
	/**
	 * The short version of the identifier of the generator
	 */
	public final static String SHORT_ID = "OPTIMISTS";
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Probability distribution of the initial/source state variables
	 */
	private NonParametric sourceStateDist;
	
	/**
	 * Remaining root samples of the {@link #sourceStateDist} to be offered to the optimizer
	 */
	private Queue<SolutionRoot> rootSamples;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param sourceStateDist	{@link #sourceStateDist}
	 * @param rootSamples		{@link #rootSamples}
	 */
	public ParticleGenerator(NonParametric sourceStateDist, Queue<SolutionRoot> rootSamples)
	{
		this.sourceStateDist	= sourceStateDist;
		this.rootSamples		= rootSamples;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	@Override
	public String getId()
	{
		return ID;
	}

	@Override
	public String getShortId()
	{
		return SHORT_ID;
	}

	@Override
	public String getParamSummary()
	{
		return "no parameters";
	}

	@Override
	public void addDiscVariable(DiscVar variable)
	{
		// Do nothing
	}

	@Override
	public void addContVariable(ContVar variable)
	{
		// Do nothing
	}

	@Override
	public void clearVariables()
	{
		// Do nothing
	}

	@Override
	public void setObjectives(ArrayList<Objective> objectives)
	{
		// Do nothing
	}

	@Override
	public ArrayList<SolutionRoot> generateSolutions(Population population, int number)
	{
		ArrayList<SolutionRoot> roots	= new ArrayList<>(number);
		
		// Retrieve root samples
		int inQueue						= rootSamples.size();
		if (inQueue > 0)
		{
			int toPoll					= Math.min(number, inQueue);
			for (int s = 0; s < toPoll; s++)
				roots.add(rootSamples.poll());
		}
		
		// Add random samples
		int toSample					= number - roots.size();
		if (toSample > 0)
		{
			double[][] samples			= sourceStateDist.sampleMultiple(toSample);
			for (int s = 0; s < samples.length; s++)
			{
				ArrayList<Double> vals	= Utilities.toArrayList(samples[s]);
				roots.add(new SolutionRoot(Assimilator.LABEL_RANDOM, null, vals, null));
			}
		}
		
		return roots;
	}

}
