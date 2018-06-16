/*
 * NormalKDEDistribution.java
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


import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Param;
import beast.core.State;
import beast.gss.TraceLog;
import beast.math.statistic.DiscreteStatistics;
import beast.util.HeapSort;
//import dr.stats.DiscreteStatistics;
//import dr.util.HeapSort;

import java.util.List;
import java.util.Random;

/**
 * @author Marc A. Suchard
 */
@Description("Distribution based on normal kernel density esitmators (and no transform on the input)")
public class NormalKDEDistribution extends KernelDensityEstimatorDistribution {
	public Input<TraceLog> traceLogInput = new Input<>("traceLog", "file containing trace log", Validate.REQUIRED);
	public Input<String> labelInput = new Input<>("label", "label of the column containing data in the trace file", Validate.REQUIRED);
	public Input<Function> xInput = new Input<>("x", "function/statistic to take distribution over");

    public static final int MINIMUM_GRID_SIZE = 512;

    @Override
    public void initAndValidate() {
    	this.p = xInput.get();
		this.traceLog = traceLogInput.get();
		this.label = labelInput.get();

		Double [] sample = traceLog.getTrace(label); 
        this.sample = new double[sample.length];
        for (int i = 0; i < sample.length; i++) {
            this.sample[i] = sample[i];
        }
        this.N = sample.length;
        lowerBound = Double.NEGATIVE_INFINITY; upperBound = Double.POSITIVE_INFINITY;
        processBounds(lowerBound, upperBound);
        setBandWidth(null);

    	
    	
        this.gridSize = MINIMUM_GRID_SIZE;
        if (this.gridSize > MINIMUM_GRID_SIZE) {
            this.gridSize = (int) Math.pow(2, Math.ceil(Math.log(this.gridSize) / Math.log(2.0)));
        }
        this.cut = 3.0;

        from = DiscreteStatistics.min(super.sample) - this.cut * this.bandWidth;
        to = DiscreteStatistics.max(super.sample) + this.cut * this.bandWidth;

        lo = from - 4.0 * this.bandWidth;
        up = to + 4.0 * this.bandWidth;

        densityKnown = false;    	super.initAndValidate();
    }
    public NormalKDEDistribution() {}
    
	public NormalKDEDistribution(TraceLog traceLog,
			String label,
			Function p) {		
		this(traceLog.getTrace(label), p);
		this.traceLog = traceLog;
		this.label = label;
		this.p = p;
		
		traceLogInput.setValue(traceLog, this);
		labelInput.setValue(label, this);
		xInput.setValue(p, this);
	}

	public NormalKDEDistribution(Double[] sample, Function p) {
        this(sample, null, null, null, p);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth, Function p) {
        this(sample, lowerBound, upperBound, bandWidth, 3.0, MINIMUM_GRID_SIZE, p);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth,
                                 int n, Function p) {
        this(sample, lowerBound, upperBound, bandWidth, 3.0, n, p);
    }

    public NormalKDEDistribution(Double[] sample, Double lowerBound, Double upperBound, Double bandWidth,
                                 double cut, int n, Function p) {
        super(sample, lowerBound, upperBound, bandWidth, p);
        this.gridSize = Math.max(n, MINIMUM_GRID_SIZE);
        if (this.gridSize > MINIMUM_GRID_SIZE) {
            this.gridSize = (int) Math.pow(2, Math.ceil(Math.log(this.gridSize) / Math.log(2.0)));
        }
        this.cut = cut;

        from = DiscreteStatistics.min(super.sample) - this.cut * this.bandWidth;
        to = DiscreteStatistics.max(super.sample) + this.cut * this.bandWidth;

        lo = from - 4.0 * this.bandWidth;
        up = to + 4.0 * this.bandWidth;

        densityKnown = false;
    }

    public double getFromPoint() {
        return from;
    }

    public double getToPoint() {
        return to;
    }

    /**
     * Returns a linear approximation evaluated at pt
     * @param x data (assumed sorted increasingly
     * @param y data
     * @param pt evaluation point
     * @param low return value if pt < x
     * @param high return value if pt > x
     * @return  evaluated coordinate
     */
    private double linearApproximate(double[] x, double[] y, double pt, double low, double high) {

        int i = 0;
        int j = x.length - 1;

        if (pt < x[i]) {
            return low;
        }
        if (pt > x[j]) {
            return high;
        }

        // Bisection search
        while (i < j - 1) {
            int ij = (i + j) / 2;
            if (pt < x[ij]) {
                j = ij;
            } else {
                i = ij;
            }
        }

        if (pt == x[j]) {
            return y[j];
        }
        if (pt == x[i]) {
            return y[i];
        }
        return y[i] + (y[j] - y[i]) * ((pt - x[i]) / (x[j] - x[i]));
    }

    private double[] rescaleAndTrim(double[] x) {
        final int length = x.length / 2;
        final double scale = 1.0 / x.length;
        double[] out = new double[length];
        for (int i = 0; i < length; ++i) {
            out[i] = x[i] * scale;
            if (out[i] < 0) {
                out[i] = 0;
            }
        }
        return out;
    }

    private double[] massdist(double[] x,
//                              double[] xmass,
                              double xlow, double xhigh, int ny) {

        int nx = x.length;
        double[] y = new double[ny * 2];

        final int ixmin = 0;
        final int ixmax = ny - 2;
        final double xdelta = (xhigh - xlow) / (ny - 1);

        for (int i = 0; i < ny; ++i) {
            y[i] = 0.0;
        }

        final double xmi = 1.0 / nx;
        for (int i = 0; i < nx; ++i) {
            final double xpos = (x[i] - xlow) /  xdelta;
            final int ix = (int) Math.floor(xpos);
            final double fx = xpos - ix;
//            final double xmi = xmass[i];

            if (ixmin <= ix && ix <= ixmax) {
                y[ix] += (1 - fx) * xmi;
                y[ix + 1] += fx * xmi;
            } else if (ix == -1) {
                y[0] += fx * xmi;
            } else if (ix == ixmax + 1) {
                y[ix] += (1 - fx) * xmi;
            }
        }
        return y;
    }

    /**
     * Override for different kernels
     * @param ordinates the points in complex space
     * @param bandWidth predetermined bandwidth
     */
    protected void fillKernelOrdinates(ComplexArray ordinates, double bandWidth) {
                final int length = ordinates.length;
        final double a = 1.0 / (Math.sqrt(2.0 * Math.PI) * bandWidth);
        final double precision = -0.5 / (bandWidth * bandWidth);
        for (int i = 0; i < length; i++) {
            final double x = ordinates.real[i];
            ordinates.real[i] = a * Math.exp(x * x * precision);
        }
    }

    protected void computeDensity() {
        makeOrdinates();
        transformData();
        densityKnown = true;
    }

    private void transformData() {
        ComplexArray Y  = new ComplexArray(massdist(this.sample, lo, up, this.gridSize));
        FastFourierTransform.fft(Y, false);

        ComplexArray product = Y.product(kOrdinates);
        FastFourierTransform.fft(product, true);

        densityPoints = rescaleAndTrim(product.real);        
    }

    private void makeOrdinates() {

        final int length = 2 * gridSize;
        if (kOrdinates == null) {
            kOrdinates = new ComplexArray(new double[length]);
        }

        // Fill with grid values
        final double max = 2.0 * (up - lo);
        double value = 0;
        final double inc = max / (length - 1);
        for (int i = 0; i <= gridSize; i++) {
            kOrdinates.real[i] = value;
            value += inc;
        }
        for (int i = gridSize + 1; i < length; i++) {
            kOrdinates.real[i] = -kOrdinates.real[length - i];
        }
        fillKernelOrdinates(kOrdinates, bandWidth);

        FastFourierTransform.fft(kOrdinates, false);
        kOrdinates.conjugate();

        // Make x grid
        xPoints = new double[gridSize];
        double x = lo;
        double delta = (up - lo) / (gridSize - 1);
        for (int i = 0; i < gridSize; i++) {
            xPoints[i] = x;
            x += delta;
        }
    }

    @Override
    protected double evaluateKernel(double x) {        
        if (!densityKnown) {
           computeDensity();
        }
        return linearApproximate(xPoints, densityPoints, x, 0.0, 0.0);
    }

    @Override
    protected void processBounds(Double lowerBound, Double upperBound) {
        if ((lowerBound != null && lowerBound != Double.NEGATIVE_INFINITY) ||
                (upperBound != null && upperBound != Double.POSITIVE_INFINITY)) {
            throw new RuntimeException("NormalKDEDistribution must be unbounded");
        }
    }

    @Override
    protected void setBandWidth(Double bandWidth) {
        if (bandWidth == null) {
            // Default bandwidth
            this.bandWidth = bandwidthNRD(sample);
        } else
            this.bandWidth = bandWidth;
                    
        densityKnown = false;
    }

//   bandwidth.nrd =
//   function (x)
//   {
//       r <- quantile(x, c(0.25, 0.75))
//       h <- (r[2] - r[1])/1.34
//       4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
//   }

    public double bandwidthNRD(double[] x) {
        int[] indices = new int[x.length];
        HeapSort.sort(x, indices);

        final double h =
                (DiscreteStatistics.quantile(0.75, x, indices) - DiscreteStatistics.quantile(0.25, x, indices)) / 1.34;
        return 1.06 *
                Math.min(Math.sqrt(DiscreteStatistics.variance(x)), h) *
                Math.pow(x.length, -0.2);
    }

    private ComplexArray kOrdinates;
    private double[] xPoints;
    private double[] densityPoints;

    private int gridSize;
    private double cut;
    private double from;
    private double to;
    private double lo;
    private double up;

    private boolean densityKnown = false;

    public static void main(String[] args) {

        long start = System.currentTimeMillis();

        Random random = new Random(1234);

        Double[] samples = new Double[10000000];
        for (int i = 0; i < samples.length; i++) {
            samples[i] = random.nextDouble();
        }
        NormalKDEDistribution nKDE = new NormalKDEDistribution(samples, null);

        for (int i = 0; i < 100; i++) {
            nKDE.evaluateKernel(random.nextDouble());
        }

        long end = System.currentTimeMillis();

        System.out.println("Time: " + (end-start));

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