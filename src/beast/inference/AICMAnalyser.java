/*
 * MarginalLikelihoodAnalysis.java
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
package beast.inference;

import beast.app.util.Application;
import beast.app.util.LogFile;
import beast.core.Description;
import beast.core.Input;
import beast.core.util.Log;
import beast.util.LogAnalyser;
import beast.util.Randomizer;

/**
 * @author Marc Suchard
 * @author Alexei Drummond
 *         <p/>
 *         Source translated from model_P.c (a component of BAli-Phy by Benjamin Redelings and Marc Suchard
 *         ported from BEAST1 MarginalLikelihoodAnalysis by Remco Bouckaert
 */

@Description("Calculate AICM of a trace log")
public class AICMAnalyser extends beast.core.Runnable {
	final public Input<LogFile> traceFileInput = new Input<>("logFile","input file containing trace log");
	public Input<Integer> burninInput = new Input<>("burnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);
	public Input<Integer> bootstrapLengthInput = new Input<>("bootstrapLength", "number of bootstrap samples used to calculate variance in AICM estimate" , 1000);
	final public Input<String> traceInput = new Input<>("trace","name of the trace to use (default likelihood)", "likelihood");

    private String traceName;
    private Double[] sample;
    //private int burnin;
    //private String analysisType; // "harmonic" for harmonic mean, "smoothed" for smoothed harmonic mean, "aicm" for AICM, "arithmetic" for arithmetic mean 
    private int bootstrapLength;

    private boolean marginalLikelihoodCalculated = false;
    private double logMarginalLikelihood;
    private double bootstrappedSE;


    
    
    @Override
    public void initAndValidate() {
    }
    
    @Override
    public void run() throws Exception {
    	LogAnalyser tracer = new LogAnalyser(traceFileInput.get().getPath(), burninInput.get(), false, false);
    	traceName = traceInput.get();
    	sample = tracer.getTrace(traceName);
    	bootstrapLength = bootstrapLengthInput.get();
    	Log.warning.println("Calculating AICM");
    	Log.warning.println("|---------|---------|---------|---------|---------|---------|---------|---------|");
    	calculate();
    	Log.info.println(toString());
    }
    
//    public MarginalLikelihoodAnalysis(double[] sample, String traceName, int burnin) {
//        this(sample, traceName, burnin, false, 1000);
//    }

    public String getTraceName() {
        return traceName;
    }

    public int getBurnin() {
        return burninInput.get();
    }

    /**
     * Constructor
     *
     * @param sample
     * @param traceName       used for 'toString' display purposes only
     * @param burnin          used for 'toString' display purposes only
     * @param analysisType
     * @param bootstrapLength a value of zero will turn off bootstrapping
     */
    public AICMAnalyser(Double[] sample, String traceName/*, int burnin, String analysisType*/, int bootstrapLength) {
        this.sample = sample;
        this.traceName = traceName;
        //this.burnin = burnin;
        //this.analysisType = analysisType;
        this.bootstrapLength = bootstrapLength;
//        System.err.println("setting burnin to "+burnin);
    }
    public AICMAnalyser() {}

    public double calculateLogMarginalLikelihood(Double[] sample) {
         return logMarginalLikelihoodAICM(sample);
    }
    
    /**
     * Calculates the AICM of a model using method-of-moments from Raftery et al. (2007)
     *
     * @param v a posterior sample of logLikelihoods
     * @return the AICM (lower values are better)
     */

    public double logMarginalLikelihoodAICM(Double[] v) {

        double sum = 0;
        final int size = v.length;
        for (int i = 0; i < size; i++)
            sum += v[i];

        double mean = sum / (double) size;

        double var = 0;
        for (int i = 0; i < size; i++)
            var += (v[i]-mean) * (v[i]-mean);
        var /= (double) size - 1;

        return 2 * var - 2 * mean;

    }

    public void calculate() {

        logMarginalLikelihood = calculateLogMarginalLikelihood(sample);
        if (bootstrapLength > 1) {
            final int sampleLength = sample.length;
            Double[] bsSample = new Double[sampleLength];
            Double[] bootstrappedLogML = new Double[bootstrapLength];
            double sum = 0;

            double progress = 0.0;
            double delta = 1.0 / bootstrapLength;
            
            //System.err.println("HME = " + logMarginalLikelihood);

            for (int i = 0; i < bootstrapLength; i++) {
            //	if (i % 10 == 0) {
            //		System.err.println((i+1) + "/" + bootstrapLength);
            //	}
            	
                fireProgress(progress);
                progress += delta;

                int[] indices = Randomizer.sampleIndicesWithReplacement(sampleLength);
                for (int k = 0; k < sampleLength; k++) {
                    bsSample[k] = sample[indices[k]];
                }
                bootstrappedLogML[i] = calculateLogMarginalLikelihood(bsSample);
                sum += bootstrappedLogML[i];
            }
            sum /= bootstrapLength;
            double bootstrappedAverage = sum;
            // Summarize bootstrappedLogML
            double var = 0;
            for (int i = 0; i < bootstrapLength; i++) {
                var += (bootstrappedLogML[i] - bootstrappedAverage) *
                        (bootstrappedLogML[i] - bootstrappedAverage);
            }
            var /= (bootstrapLength - 1.0);
            bootstrappedSE = Math.sqrt(var);
        }

        fireProgress(1.0);

        marginalLikelihoodCalculated = true;
    }



    public double getLogMarginalLikelihood() {
        if (!marginalLikelihoodCalculated) {
            calculate();
        }
        return logMarginalLikelihood;
    }

    public double getBootstrappedSE() {
        if (!marginalLikelihoodCalculated) {
            calculate();
        }
        return bootstrappedSE;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        
        //if (analysisType.equals("smoothed")) {
        //    sb.append("log marginal likelihood (using smoothed harmonic mean)");
        //}
        //else if (analysisType.equals("aicm")) {
            sb.append("AICM");
        //}
        //else if (analysisType.equals("arithmetic")) {
        //	sb.append("log marginal likelihood (using arithmetic mean)");
        //}
        //else {
        //    sb.append("log marginal likelihood (using harmonic mean)");
        //}
        
        sb.append(" from ")
        		.append(traceFileInput.get().getName())
        		.append(":")
                .append(traceName)
                .append(" = ")
                .append(String.format("%5.4f", getLogMarginalLikelihood()));
        if (bootstrapLength > 1) {
            sb.append(" +/- ")
                    .append(String.format("%5.4f", getBootstrappedSE()));
        } else {
            sb.append("           ");
        }

        sb.append(" burnin=").append(burninInput.get());
        if (bootstrapLength > 1)
            sb.append(" replicates=").append(bootstrapLength);
//        sb.append("\n");

        return sb.toString();

    }
//
//
//    private TaskListener listener = null;
//
//    public void setTaskListener(TaskListener listener) {
//        this.listener = listener;
//    }
//
    int done = 0;
    private void fireProgress(double progress) {
    	while (done < progress * 80) {
    		Log.warning.print("*");
    		done++;
    	}
    	if (progress == 1.0) {
    		Log.warning.println();
    	}
    }

    
    public static void main(String[] args) throws Exception {
		new Application(new AICMAnalyser(), "AICM analyser", args);
	}
}
