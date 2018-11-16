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

import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Hashtable;

import optimists.Particle;
import vic.Soil;
import vic.routing.MuskingumNetwork;
import vic.routing.State;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class VICParticle implements Particle
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	private VICAssimilator				parent;
	private String						id;
	private ArrayList<Double>			initialArray;
	private ArrayList<Double>			finalArray;
	private State						initialState;
	private State						finalState;
	private ArrayList<Soil>				soils;
	private MuskingumNetwork			channels;
	private Hashtable<String, Double>	directFractions;
	private ArrayList<Double>			streamflow;
	private ArrayList<Double>			fitness;
	private boolean						valid;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------

	public VICParticle(VICAssimilator parent, String id, ArrayList<Double> initialArray, 
			ArrayList<Double> finalArray, State initialState, State finalState, 
			ArrayList<Soil> soils, MuskingumNetwork channels, 
			Hashtable<String, Double> directFractions, ArrayList<Double> streamflow,
			ArrayList<Double> fitness, boolean valid)
	{
		this.parent				= parent;
		this.id					= id;
		this.initialArray		= initialArray;
		this.finalArray			= finalArray;
		this.initialState		= initialState;
		this.finalState			= finalState;
		this.soils				= soils;
		this.channels			= channels;
		this.directFractions	= directFractions;
		this.streamflow			= streamflow;
		this.fitness			= fitness;
		this.valid				= valid;
	}

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------

	@Override
	public Particle createNew(int index, LocalDateTime start, LocalDateTime end,
								ArrayList<Double> initialStateArray)
	{
		return parent.createNewParticle(index, start, end, initialStateArray);
	}

	@Override
	public String getId()
	{
		return id;
	}

	@Override
	public ArrayList<Double> getSourceStateArray()
	{
		return initialArray;
	}

	@Override
	public ArrayList<Double> getTargetStateArray()
	{
		return finalArray;
	}

	public State getInitialState()
	{
		return initialState;
	}

	public State getFinalState()
	{
		return finalState;
	}

	public ArrayList<Soil> getSoils()
	{
		return soils;
	}

	public MuskingumNetwork getChannels()
	{
		return channels;
	}

	public Hashtable<String, Double> getDirectFractions()
	{
		return directFractions;
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

	@Override
	public boolean isValid()
	{
		return valid;
	}

}
