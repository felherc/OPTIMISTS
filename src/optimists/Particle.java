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

/**
 * Interfaces the data assimilation algorithm with the user's modeling framework by defining the
 * methods that a candidate configuration of the model should expose. The implementing class is 
 * responsible for evaluating candidate configurations of the model and for sharing the performance 
 * metrics of that configuration. 
 * <br><br>Such implementation or "particle" is expected to include an initial "source" state that 
 * is representable by an array of continuous values (double-precision), a similar final "target" 
 * state after an arbitrary time step, and a series of performance metrics that allow comparing it 
 * to other particles (given the user-selected objectives). The class implementing the interface is 
 * responsible for computing the target state and the performance metrics based on a provided 
 * source state. 
 * <br><br>These particles are used in the context of the data assimilation framework to allow 
 * representing multivariate probability distributions for the state variables using kernel density 
 * estimation.
 * 
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public interface Particle
{

	// --------------------------------------------------------------------------------------------
	// Methods
	// --------------------------------------------------------------------------------------------
	
	/**
	 * Creates a new particle object with the provided source state, runs the model, and computes 
	 * both the target state and the fitness variables. The provided values must be stored to be 
	 * able to return them using the method {@link #getSourceStateArray()}. These values may be 
	 * modified if necessary (when special constraints must be met, for instance). The returned
	 * object must be able to return valid results for the following methods: {@link #getId()}, 
	 * {@link #getSourceStateArray()}, {@link #getTargetStateArray()}, {@link #getReport()},
	 * {@link #getFitness()}, and {@link #compareTo(int, Particle)}.
	 * @param index				A consecutive integer value to be optionally used in the 
	 * 							identification of the new solution. The identification is returned 
	 * 							by the {@link #getId()} method.
	 * @param start				Date and time corresponding to the provided initial or source state
	 * @param end				Date and time corresponding to the requested final or target state
	 * @param sourceStateArray	Array of values that represent the source state. The values are
	 * 							organized according to the list provided at 
	 * 					{@link OPTIMISTS#ParetoDA(String, int, Particle, ArrayList, ArrayList, int)}
	 * @return The new particle object
	 */
	public Particle createNew(int index, LocalDateTime start, LocalDateTime end, 
								ArrayList<Double> sourceStateArray);
	
	/**
	 * Returns the identifier of the particle. Must not be an empty string (""). Each particle must
	 * have a unique identifier.
	 * @return The identifier of the particle
	 */
	public String getId();
	
	/**
	 * @return Array of values that represent the initial or source state
	 */
	public ArrayList<Double> getSourceStateArray();
	
	/**
	 * @return Array of values that represent the final or target state
	 */
	public ArrayList<Double> getTargetStateArray();
	
	/**
	 * Returns a string with the header of the columns in the report table. The report should
	 * include the fitness values for the multiple optimization objectives and can additionally 
	 * contain any other values of interest. The individual header strings for each field should be 
	 * separated by a tab (\t) character.
	 * @return A string with the header of the columns in the report table
	 */
	public String getReportHeader();
	
	/**
	 * Returns the values of the particle to be printed in the report table. The report should
	 * include the fitness values for the multiple optimization objectives and can additionally 
	 * contain any other values of interest. The individual strings for each field should be 
	 * separated by a tab (\t) character. The fields should correspond to those in the header 
	 * string returned by the {@link #getReportHeader()} method.
	 * @return The values of the particle to be printed in the report table
	 */
	public String getReport();
	
	/**
	 * @param objective Index of the optimization objective
	 * @return The fitness value for the requested optimization objective
	 */
	public double getFitness(int objective);
	
	/**
	 * Compares this particle with another one according to a specific optimization objective. 
	 * Returns a negative integer, zero, or a positive integer as this solution is less fitted, 
	 * equally fitted, or more fitted than the other according to the provided objective. Will only 
	 * be called for comparisons with objectives that do not use a single fitness value (as 
	 * provided by the {@link #getFitness} method).
	 * @param objective Index of the optimization objective
	 * @param other The other solution to compare
	 * @return A negative integer, zero, or a positive integer as this solution is less fitted, 
	 * equally fitted, or more fitted than the other
	 */
	public int compareTo(int objective, Particle other);

	/**
	 * @return <code>true</code> if the particle was successfully created and if it should be
	 * considered for the ensemble. <code>false</code> if for some reason the particle could not be 
	 * successfully created and should not be considered for the ensemble. In the latter case, the
	 * particle will not count towards the candidate count limit.
	 */
	public boolean isValid();
	
}