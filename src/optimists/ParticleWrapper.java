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

import java.time.LocalDateTime;
import java.util.ArrayList;

import maestro_mo.solution.Solution;
import probDist.multiVar.tools.ContMultiSample;

/**
 * Wraps a {@link Particle} implementation to include additional information such as the index, the 
 * start and end time stamps of the simulation, and the weight.
 * 
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class ParticleWrapper implements Particle, Solution, ContMultiSample
{
	
	// --------------------------------------------------------------------------------------------
	// Attributes
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Unique index of the particle. -1 if not index has been assigned yet.
	 */
	private int index;
	
	/**
	 * The particle being wrapped
	 */
	private Particle particle;
	
	/**
	 * Time stamp corresponding to the source state of the particle
	 */
	private LocalDateTime start;
	
	/**
	 * Time stamp corresponding to the target state of the particle
	 */
	private LocalDateTime end;
	
	/**
	 * The relative weight of the particle used to define a multivariate probability distribution
	 * of the target state using kernel density estimation
	 */
	private double weight;
	
	// --------------------------------------------------------------------------------------------
	// Constructor
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @param particle		{@link #particle}
	 * @param start			{@link #start}
	 * @param end			{@link #end}
	 */
	public ParticleWrapper(Particle particle, LocalDateTime start, LocalDateTime end)
	{
		this.index		= -1;
		this.particle	= particle;
		this.start		= start;
		this.end		= end;
		this.weight		= Double.NaN;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * @return {@link #index}
	 */
	public int getIndex()
	{
		return index;
	}

	/**
	 * @param index {@link #index}
	 */
	public void setIndex(int index)
	{
		this.index = index;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods - from Particle
	// --------------------------------------------------------------------------------------------

	@Override
	public Particle createNew(int index, LocalDateTime start, LocalDateTime end, 
								ArrayList<Double> sourceStateArray)
	{
		return particle.createNew(index, start, end, sourceStateArray);
	}

	@Override
	public String getId()
	{
		return particle.getId();
	}
	
	/**
	 * @return {@link #particle}
	 */
	public Particle getParticle()
	{
		return particle;
	}

	/**
	 * @return {@link #start}
	 */
	public LocalDateTime getStart()
	{
		return start;
	}

	/**
	 * @return {@link #end}
	 */
	public LocalDateTime getEnd()
	{
		return end;
	}

	@Override
	public ArrayList<Double> getSourceStateArray()
	{
		return particle.getSourceStateArray();
	}

	@Override
	public ArrayList<Double> getTargetStateArray()
	{
		return particle.getTargetStateArray();
	}

	@Override
	public String getReportHeader()
	{
		return particle.getReportHeader();
	}

	@Override
	public String getReport()
	{
		return particle.getReport();
	}

	@Override
	public double getFitness(int objective)
	{
		return particle.getFitness(objective);
	}

	@Override
	public int compareTo(int objective, Particle other)
	{
		return particle.compareTo(objective, other);
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods - from ContMultiSample
	// --------------------------------------------------------------------------------------------

	@Override
	public double getWeight()
	{
		return weight;
	}

	@Override
	public void setWeight(double weight)
	{
		this.weight = weight;
	}

	@Override
	public ArrayList<Double> getValues()
	{
		return particle.getTargetStateArray();
	}

	@Override
	public ContMultiSample copy()
	{
		// TODO Implement
		return null;
	}
	
	// --------------------------------------------------------------------------------------------
	// Methods - from Solution
	// --------------------------------------------------------------------------------------------

	@Override
	public Solution createNew(int index, ArrayList<Integer> discValues, 
								ArrayList<Double> contValues, Object extra)
	{
		return new ParticleWrapper(particle.createNew(index, start, end, contValues), start, end);
	}

	@Override
	public ArrayList<Integer> getDiscValues()
	{
		return null;
	}

	@Override
	public ArrayList<Double> getContValues()
	{
		return particle.getSourceStateArray();
	}

	@Override
	public boolean isValid()
	{
		if (particle == null)
			return false;
		return particle.isValid();
	}

	@Override
	public int compareTo(int objective, Solution other)
	{
		return particle.compareTo(objective, (Particle) other);
	}

	@Override
	public boolean optimizationConverged()
	{
		return false;
	}
	
}
