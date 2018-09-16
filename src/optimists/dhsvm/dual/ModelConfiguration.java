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

package optimists.dhsvm.dual;

import java.util.ArrayList;

import dhsvm.Soil;
import dhsvm.Vegetation;
import dhsvm.grid.State;
import dhsvm.stream.StreamNetwork;

/**
 * Citation: Hernández, F., & Liang, X. (2018). "Hybridizing Bayesian and variational data
 * assimilation for high-resolution hydrologic forecasting". Hydrol. Earth Syst. Sci.
 * (Please cite this article if you use OPTIMISTS.)
 * 
 * @author Felipe Hernández (developer)
 * @author Xu Liang (advisor)
 */
public class ModelConfiguration
{

	public ArrayList<Soil>			soils;
	public ArrayList<Vegetation>	vegetations;
	public StreamNetwork			network;
	public State					initialState;
	
	public ModelConfiguration(ArrayList<Soil> soils, ArrayList<Vegetation> vegetations,
								StreamNetwork network, State initialState)
	{
		this.soils			= soils;
		this.vegetations	= vegetations;
		this.network		= network;
		this.initialState	= initialState;
	}
	
}
