/*
 * MultivariateKDEDistribution.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package beast.gss.distributions;


import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Param;
import beast.core.State;

/**
 * @author Guy Baele
 */
@Description("Multivariate kernel density esitmators that assumes input variable has independent components")
public class MultivariateKDEDistribution extends Distribution {
	public Input<List<KernelDensityEstimatorDistribution>> distInput = new Input<>("dist","distributions, one for each of the dimensions of the parameter", new ArrayList<>());
	public Input<Function> xInput = new Input<>("x", "parameter to which this distribution applies");
	
	public static final String TYPE = "multivariateKDE";
    public static final boolean DEBUG = false;
	
	private KernelDensityEstimatorDistribution[] multivariateKDE;
	private Function p;
	private int dimension;
	//private boolean[] flags;
	

	@Override
	public void initAndValidate() {
		this.p = xInput.get();
		this.multivariateKDE = distInput.get().toArray(new KernelDensityEstimatorDistribution[]{});
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE)");
		}
		this.dimension = multivariateKDE.length;
		super.initAndValidate();
	}
	
	public MultivariateKDEDistribution () {}
	
	public MultivariateKDEDistribution (KernelDensityEstimatorDistribution[] multivariateKDE, Function p) {
		this.p = p;
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		
		xInput.setValue(p, this);
		for (KernelDensityEstimatorDistribution k : multivariateKDE) {
			distInput.get().add(k);
		}
		/*for (int i = 0; i < dimension; i++) {
			flags[i] = true;
		}*/

	}
	
	public MultivariateKDEDistribution (KernelDensityEstimatorDistribution[] multivariateKDE, boolean[] flags) {
		
		if (multivariateKDE.length <= 0) {
			throw new RuntimeException("Creation error in MultivariateKDEDistribution(Distribution[] multivariateKDE, boolean[] flags)");
		}
		
		this.multivariateKDE = multivariateKDE;
		this.dimension = multivariateKDE.length;
		//this.flags = flags;

		xInput.setValue(p, this);
		for (KernelDensityEstimatorDistribution k : multivariateKDE) {
			distInput.get().add(k);
		}
	}

	@Override
	public double calculateLogP() {
		logP = 0;
		for (int i = 0; i < p.getDimension(); i++) {
			logP += multivariateKDE[i].logPdf((double) p.getArrayValue(i));
		}
		return logP;
	}
	
	public double logPdf(double[] x) {
		
		double logPdf = 0;
		
		if (x.length != dimension) {
            throw new IllegalArgumentException("data array is of the wrong dimension");
        }
		
		for (int i = 0; i < dimension; i++) {
			//if (flags[i]) {
			logPdf += multivariateKDE[i].logPdf(x[i]);
			//}
		}

        if (DEBUG){
            System.err.println("MultivariateKDEDistribution, dimension = " + dimension);
            for (int i = 0; i < dimension; i++) {
                System.err.println(i + ", " + "x[i] = " + x[i] + ", logPdf = " + multivariateKDE[i].logPdf(x[i]));
                //System.err.println("    mean = " + multivariateKDE[i].mean() + ", variance = " + multivariateKDE[i].variance());
            }
        }
		
		return logPdf;
	}

	public double[][] getScaleMatrix() {
		throw new RuntimeException("Not yet implemented");
	}

	public double[] getMean() {
		throw new RuntimeException("Not yet implemented");
	}

	public String getType() {
		return TYPE;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
}
