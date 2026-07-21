/*
 * KernelDensityEstimatorDistribution.java
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

package modelselection.gss.distribution;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.parameter.Parameter;
import modelselection.gss.TraceLog;

//import dr.math.UnivariateFunction;

/**
 * @author Marc Suchard
 */
@Description("Distribution based on kernel density esitmators")
public abstract class KernelDensityEstimatorDistribution extends Distribution {	
	protected Function p;
	protected TraceLog traceLog;
	protected String label;

    public KernelDensityEstimatorDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth, Function p) {
    	this.p = p;
        this.sample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.sample[i] = sample[i];
        }
        this.N = sample.length;
        processBounds(lowerBound, upperBound);
        setBandWidth(bandWidth);
    }

    public KernelDensityEstimatorDistribution() {
		// TODO Auto-generated constructor stub
	}

	abstract protected double evaluateKernel(double x);

    abstract protected void processBounds(Double lowerBound, Double upperBound);

    abstract protected void setBandWidth(Double bandWidth);

    /**
     * probability density function of the distribution
     *
     * @param x argument
     * @return pdf value
     */
    public double pdf(double x) {
        return evaluateKernel(x);
    }

    /**
     * the natural log of the probability density function of the distribution
     *
     * @param x argument
     * @return log pdf value
     */
    public double logPdf(double x) {
        return Math.log(pdf(x));
    }
    
	@Override
	public double calculateLogP() {
		logP = 0;
		if (p instanceof Parameter) {
			for (int i = 0; i < p.getDimension(); i++) {
				if (p.getArrayValue(i) < (double) ((Parameter)p).getLower()) {
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}
				if (p.getArrayValue(i) > (double) ((Parameter)p).getUpper()) {
					logP = Double.NEGATIVE_INFINITY;
					return logP;
				}			
			}
		}

		for (int i = 0; i < p.getDimension(); i++) {
			logP += logPdf(p.getArrayValue(i));
		}
		return logP;
	}

    /**
     * cumulative density function of the distribution
     *
     * @param x argument
     * @return cdf value
     */
    public double cdf(double x) {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * quantile (inverse cumulative density function) of the distribution
     *
     * @param y argument
     * @return icdf value
     */
    public double quantile(double y) {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * mean of the distribution
     *
     * @return mean
     */
    public double mean() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * variance of the distribution
     *
     * @return variance
     */
    public double variance() {
        throw new RuntimeException("Not Implemented.");
    }

    /**
     * @return a probability density function representing this distribution
     */
//    public UnivariateFunction getProbabilityDensityFunction() {
//        throw new RuntimeException("Not Implemented.");
//    }

    public double getBandWidth() {
        return bandWidth;
    }

    public enum Type {
        GAUSSIAN("Gaussian"),
        GAMMA("Gamma"),
        LOGTRANSFORMEDGAUSSIAN("LogTransformedGaussian"),
        BETA("Beta");

        private Type(String text) {
            this.text = text;
        }

        public String getText() {
            return text;
        }

        public static Type parseFromString(String text) {
            for (Type format : Type.values()) {
                if (format.getText().compareToIgnoreCase(text) == 0)
                    return format;
            }
            return null;
        }

        private final String text;
    }

    protected int N;
    protected double lowerBound;
    protected double upperBound;
    protected double bandWidth;
    protected double[] sample;
    
    @Override
    protected boolean requiresRecalculation() {
    	return true;
    }


	public TraceLog getTraceLog() {
		return traceLog;
	}
	public void setTraceLog(TraceLog traceLog) {
		this.traceLog = traceLog;
	}
	public String getLabel() {
		return label;
	}
	public void setLabel(String label) {
		this.label = label;
	}
	public Function getX() {
		return p;
	}
	public void setX(Function p) {
		this.p = p;
	}
}
