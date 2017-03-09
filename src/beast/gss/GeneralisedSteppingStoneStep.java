package beast.gss;




import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.app.util.LogFile;
import beast.app.util.TreeFile;
import beast.core.*;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.cpo.BEASTRunAnalyser;
import beast.evolution.tree.TreeDistribution;
import beast.gss.distributions.GSSTreeDistribution;
import beast.gss.distributions.KernelDensityEstimatorDistribution;
import beast.gss.distributions.MultivariateKDEDistribution;
import beast.gss.distributions.NormalKDEDistribution;
import beast.inference.PathSamplingStep;
import beast.math.distributions.MRCAPrior;
import beast.util.LogAnalyser;
import beast.util.Randomizer;

@Description("Calculate marginal likelihood through generalised path sampling for a single step")
public class GeneralisedSteppingStoneStep extends PathSamplingStep {
	final public Input<LogFile> traceFileInput = new Input<>("logFile","input file containing trace log. If not specified, use trace log file from XML");
	final public Input<TreeFile> treeFileInput = new Input<>("treeFile","input file containing tree log. If not specified, use tree log file from XML");
	public Input<Integer> traceBurninInput = new Input<>("traceBurnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);

	@Override
	public void initAndValidate() {

		LogAnalyser tracelog;
		try {
			tracelog = BEASTRunAnalyser.getTraceLog(this, traceFileInput.get(), traceBurninInput.get());
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
		CompoundDistribution prior = getPrior(this);
		CompoundDistribution altPrior = getAltPrior(prior, tracelog);

		Distribution posterior = posteriorInput.get();
		if (posterior instanceof CompoundDistribution) {
			CompoundDistribution p = (CompoundDistribution) posterior;
			p.pDistributions.get().add(altPrior);
		} else {
			throw new IllegalArgumentException("Expected CompoundDistribution for distribution input");
		}

		super.initAndValidate();
	}
	
	private CompoundDistribution getAltPrior(CompoundDistribution prior, LogAnalyser tracelog) {
		List<Distribution> altPrior = new ArrayList<>();
		for (Distribution d : prior.pDistributions.get()) {
			if (d instanceof TreeDistribution) {
				Distribution altTreeDist = getAltTreeDist((TreeDistribution) d);
				altPrior.add(altTreeDist);
			} else if (d instanceof MRCAPrior) {
				altPrior.add(d);
			} else if (d instanceof beast.math.distributions.Prior) {
				Distribution altPriorDist = getAltPriorDist((beast.math.distributions.Prior) d, tracelog);
				altPrior.add(altPriorDist);
			} else {
				throw new IllegalArgumentException("Don't know how to handle distributio " + d.getID() + " of type " + d.getClass().getName());
			}
			
		}
		CompoundDistribution cd = new CompoundDistribution();
		cd.initByName("distribution", altPrior);
		cd.setID("GSSPrior");
		return cd;
	}
	
	private Distribution getAltTreeDist(TreeDistribution d) {
		GSSTreeDistribution ccDistr = new GSSTreeDistribution(treeFileInput.get(), d.treeInput.get(), traceBurninInput.get());
		CompoundDistribution cd = new CompoundDistribution();
		cd.initByName("distribution", ccDistr);
		return cd;
	}

	private Distribution getAltPriorDist(beast.math.distributions.Prior d, LogAnalyser tracelog) {
    	List<String> labels = tracelog.getLabels();
    	Double[][] trace = null;
		Function f = d.m_x.get();
		String id = ((BEASTInterface)f).getID();
		String shortid = id.contains(".") ? id.substring(0, id.lastIndexOf('.')): id;
		id += ".1";
		int index = labels.indexOf(shortid);
		if (index < 0) {
			shortid += ".1";
			index = labels.indexOf(shortid);
		}
		if (index >= 0) {
			trace = new Double[f.getDimension()][];
			for (int j = 0; j < f.getDimension(); j++) {
				Double [] v = tracelog.getTrace(1 + index + j);
				trace[j] = v;
			}
		} else {
			Log.warning.println("Did not find " + id + " in log.");
		}
    	Distribution altDist = null;
    	
		int dim = f.getDimension();
    	if (dim == 1) {
    		altDist = new NormalKDEDistribution(trace[0], f); 
    	} else {
    		KernelDensityEstimatorDistribution [] multivariateKDE = new KernelDensityEstimatorDistribution[dim]; 
    		for (int i = 0; i < dim; i++) {
    			multivariateKDE[i] = new NormalKDEDistribution(trace[i], null);
    		}
    		altDist = new MultivariateKDEDistribution(multivariateKDE, f);
    	}
    			//EmpiricalDist(d.m_x.get(), trace);
		
		return altDist;
	}

	private CompoundDistribution getPrior(MCMC mcmc) {
		if (!(mcmc.posteriorInput.get() instanceof CompoundDistribution)) {
			throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with posterior CompoundDistribution");
		}
		CompoundDistribution posterior = (CompoundDistribution) mcmc.posteriorInput.get();
		for (Distribution d : posterior.pDistributions.get()) {
			if (d.getID().equals("prior")) {
				if (!(d instanceof CompoundDistribution)) {
					throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with prior CompoundDistribution");
				}
				return (CompoundDistribution) d;
			}
		}
		throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with a 'prior' in the posterior");
	} // getPrior


	 /**
     * main MCMC loop *
     */
	@Override
    protected void doLoop() {
        oldLogLikelihood = pDists[pDists.length - 1].calculateLogP() * (1.0 - beta); // GSS prior
        for (int i = 0; i < pDists.length; i++) { // posterior
            oldLogLikelihood += pDists[i].calculateLogP() * beta;
		}

        for (int iSample = -burnIn; iSample <= chainLength; iSample++) {
            final int currentState = iSample;

            state.store(currentState);
            if (storeEvery > 0 && iSample % storeEvery == 0 && iSample > 0) {
                state.storeToFile(iSample);
                // Do not store operator optimisation information
                // since this may not be valid for the next step
                // especially when sampling from the prior only
            	//operatorSchedule.storeToFile();
            }

            Operator operator = operatorSchedule.selectOperator();
            //System.out.print("\n" + iSample + " " + operator.getName()+ ":");

            final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
            Evaluator evaluator = null;

            if (evaluatorDistribution != null) {
                evaluator = new Evaluator() {
                    @Override
                    public double evaluate() {
                        double logP = 0.0;

                        state.storeCalculationNodes();
                        state.checkCalculationNodesDirtiness();

                        try {
                            logP = evaluatorDistribution.calculateLogP();
                        } catch (Exception e) {
                            e.printStackTrace();
                            System.exit(1);
                        }

                        state.restore();
                        state.store(currentState);

                        return logP;
                    }
                };
            }

            double fLogHastingsRatio = operator.proposal(evaluator);

            if (fLogHastingsRatio != Double.NEGATIVE_INFINITY) {

                state.storeCalculationNodes();
                state.checkCalculationNodesDirtiness();

                posterior.calculateLogP();
                
                newLogLikelihood = pDists[pDists.length - 1].getArrayValue() * (1.0 - beta); // GSS prior
                for (int i = 0; i < pDists.length; i++) { // posterior
                	newLogLikelihood += pDists[i].getArrayValue() * beta;
        		}

                logAlpha = newLogLikelihood - oldLogLikelihood + fLogHastingsRatio; //CHECK HASTINGS
                //System.out.println(logAlpha + " " + fNewLogLikelihood + " " + fOldLogLikelihood);
                if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                    // accept
                    oldLogLikelihood = newLogLikelihood;
                    state.acceptCalculationNodes();

                    if (iSample >= 0) {
                        operator.accept();
                    }
                    //System.out.print(" accept");
                } else {
                    // reject
                    if (iSample >= 0) {
                        operator.reject();
                    }
                    state.restore();
                    state.restoreCalculationNodes();
                    //System.out.print(" reject");
                }
                state.setEverythingDirty(false);
            } else {
                // operation failed
                if (iSample >= 0) {
                    operator.reject();
                }
                state.restore();
                //System.out.print(" direct reject");
            }
            log(iSample);

            if (iSample >= 0) {
            	operator.optimize(logAlpha);
            }
            callUserFunction(iSample);
        }
    }
}
