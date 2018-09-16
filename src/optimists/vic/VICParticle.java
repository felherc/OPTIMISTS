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

import java.time.LocalDateTime;
import java.util.ArrayList;

import optimists.Particle;
import optimists.vic.VICAssimilator;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICParticle implements Particle
{

	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private VICAssimilator		parent;
	private String				id;
	private ArrayList<Double>	sourceStateArray;
	private ArrayList<Double>	targetStateArray;
	private State				sourceState;
	private State				targetState;
	private ArrayList<Double>	streamflow;
	private ArrayList<Double>	fitness;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	public VICParticle(VICAssimilator parent, String id, ArrayList<Double> sourceStateArray,
			ArrayList<Double> targetStateArray, State sourceState, State targetState, 
			ArrayList<Double> streamflow, ArrayList<Double> fitness)
	{
		this.parent				= parent;
		this.id					= id;
		this.sourceStateArray	= sourceStateArray;
		this.targetStateArray	= targetStateArray;
		this.sourceState		= sourceState;
		this.targetState		= targetState;
		this.streamflow			= streamflow;
		this.fitness			= fitness;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	@Override
	public Particle createNew(int index, LocalDateTime start, LocalDateTime end, 
								ArrayList<Double> sourceStateArray)
	{
		return parent.createNewParticle(index, start, end, sourceStateArray);
	}

	@Override
	public String getId()
	{
		return id;
	}

	@Override
	public ArrayList<Double> getSourceStateArray()
	{
		return sourceStateArray;
	}

	@Override
	public ArrayList<Double> getTargetStateArray()
	{
		return targetStateArray;
	}
	
	public State getSourceState()
	{
		return sourceState;
	}
	
	public State getTargetState()
	{
		return targetState;
	}
	
	@Override
	public boolean isValid()
	{
		return true;
	}

	public ArrayList<Double> getStreamflow()
	{
		return streamflow;
	}

	@Override
	public String getReportHeader()
	{
		return parent.getParticleReportHeader();
	}

	@Override
	public String getReport()
	{
		return parent.getParticleReport(this);
	}

	@Override
	public double getFitness(int objective)
	{
		return fitness.get(objective);
	}

	@Override
	public int compareTo(int objective, Particle other)
	{
		return 0;
	}

}
