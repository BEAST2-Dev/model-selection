package modelselection.gss;


import beast.core.*;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.util.Randomizer;
import modelselection.inference.PathSamplingStep;

@Citation(value= "Xie, W., Lewis, P. O., Fan, Y., Kuo, L., & Chen, M. H. (2011).\n" +
        "Improving marginal likelihood estimation for Bayesian phylogenetic model selection.\n" +
        "Systematic biology, 60(2), 150-160."
        , year = 2011, firstAuthorSurname = "Xie",
        DOI="10.1093/sysbio/syq085")
@Description("Calculate marginal likelihood through generalised path sampling for a single step")
public class GeneralisedSteppingStoneStep extends PathSamplingStep {
//	final public Input<LogFile> traceFileInput = new Input<>("logFile","input file containing trace log. If not specified, use trace log file from XML");
//	final public Input<TreeFile> treeFileInput = new Input<>("treeFile","input file containing tree log. If not specified, use tree log file from XML");
//	public Input<Integer> traceBurninInput = new Input<>("traceBurnin", "percentage of the log file to disregard as burn-in (default 10)" , 10);

	final public Input<Distribution> samplingDistributionInput = new Input<>("samplingDistribution",
			"probability distribution to sample from (e.g. a prior). "
					+ "If not specified, everything but the likelihood will be used as sampling distribution.");
	
	Distribution samplingDistribution;

	@Override
	public void initAndValidate() {
		samplingDistribution = samplingDistributionInput.get();		
		if (samplingDistribution != null) {
			CompoundDistribution d = (CompoundDistribution) posteriorInput.get();
			d.pDistributions.setValue(samplingDistribution, d);
		} else {
			samplingDistribution = posteriorInput.get();
		}
		
		super.initAndValidate();
	}

	 /**
     * main MCMC loop *
     */
	@Override
    protected void doLoop() {
		robustlyCalcPosterior(posterior);
        oldLogLikelihood = pDists[pDists.length - 1].getArrayValue() * (1.0 - beta); // GSS prior
        for (int i = 0; i < pDists.length - 1; i++) { // posterior
            oldLogLikelihood += pDists[i].getArrayValue() * beta;
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
                
                newLogLikelihood = samplingDistribution.getArrayValue() * (1.0 - beta); // GSS prior
                for (int i = 0; i < pDists.length - 1; i++) { // posterior
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
