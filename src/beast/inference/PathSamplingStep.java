package beast.inference;

import java.io.IOException;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.StateNodeInitialiser;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.util.Randomizer;



@Description("Calculate marginal likelihood through path sampling for a single step")
public class PathSamplingStep extends MCMC {

	public Input<Double> betaInput = new Input<Double>("beta","power used for likelihood: 1 = using full posterior, 0 = using prior only", 1.0);

	double beta;
    Distribution prior;
    Distribution msc;
	Distribution likelihood;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		beta = betaInput.get();
		posterior = posteriorInput.get();
		// expect compound distribution with likelihood and prior
		if (!(posterior instanceof CompoundDistribution)) {
			throw new IllegalArgumentException("Expected posterior being a CompoundDistribution");
		}
		CompoundDistribution d = (CompoundDistribution) posterior;
		List<Distribution> list = d.pDistributions.get();
        if (list.size() > 3 || list.size() < 2) {
            throw new IllegalArgumentException("Expected one prior, one likelihood and an optional multispecies coalescent distribution.");
        }

        prior = null;
        msc = null;
        likelihood = null;
        for (Distribution subDist: list) {
            final String distID = subDist.getID().toLowerCase();
            if (distID.startsWith("prior")) prior = subDist;
            if (distID.startsWith("speciescoalescent")) msc = subDist;
            if (distID.startsWith("likelihood")) likelihood = subDist;
        }

        if (msc == null && list.size() == 3)
            throw new IllegalArgumentException("Expected one prior, one likelihood and an optional multispecies coalescent distribution.");

        if (prior == null)
            throw new IllegalArgumentException("No prior distribution found in posterior.");
        if (likelihood == null)
            throw new IllegalArgumentException("No likelihood distribution found in posterior.");

        loggers = loggersInput.get();
	}
	
	
    @Override
    public void run() throws SAXException, IOException, ParserConfigurationException {
        // set up state (again). Other plugins may have manipulated the
        // StateNodes, e.g. set up bounds or dimensions
        state.initAndValidate();
        // also, initialise state with the file name to store and set-up whether to resume from file
        state.setStateFileName(stateFileName);
        operatorSchedule.setStateFileName(stateFileName);

        burnIn = burnInInput.get();
        chainLength = chainLengthInput.get();
        int nInitiliasiationAttemps = 0;
        state.setEverythingDirty(true);
        posterior = posteriorInput.get();

        if (restoreFromFile) {
            state.restoreFromFile();
            operatorSchedule.restoreFromFile();
            burnIn = 0;
            oldLogLikelihood = robustlyCalcPosterior(posterior);
        } else {
            do {
                for (StateNodeInitialiser initialiser : initialisersInput.get()) {
                    initialiser.initStateNodes();
                }
                oldLogLikelihood = robustlyCalcPosterior(posterior);
            } while (Double.isInfinite(oldLogLikelihood) && nInitiliasiationAttemps++ < 10);
        }
        long tStart = System.currentTimeMillis();

        // do the sampling
        logAlpha = 0;
        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));


//        System.err.println("Start state:");
//        System.err.println(state.toString());

        System.err.println("Start likelihood: " + oldLogLikelihood + " " + (nInitiliasiationAttemps > 1 ? "after " + nInitiliasiationAttemps + " initialisation attempts" : ""));
        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
            reportLogLikelihoods(posterior, "");
            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.");
        }

        // initialises log so that log file headers are written, etc.
        for (Logger log : loggersInput.get()) {
            log.init();
        }

        doLoop();

        operatorSchedule.showOperatorRates(System.out);
        long tEnd = System.currentTimeMillis();
        System.out.println("Total calculation time: " + (tEnd - tStart) / 1000.0 + " seconds");
        close();

        System.err.println("End likelihood: " + oldLogLikelihood);
//        System.err.println(state);
        state.storeToFile(chainLength);
        // Do not store operator optimisation information
        // since this may not be valid for the next step
        // especially when sampling from the prior only
//        operatorSchedule.storeToFile();
    } // run;
	
	
	
    /**
     * main MCMC loop *
     */
    protected void doLoop() {
    	
        double logPriorProb = prior.calculateLogP();
        double logMSC = (msc == null) ? 0.0 : msc.calculateLogP();
        double logLikelihood = likelihood.calculateLogP();
        oldLogLikelihood = logPriorProb + logMSC + (logLikelihood * beta);
    	
    	
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
                logPriorProb = prior.getArrayValue();
                logMSC = (msc == null) ? 0.0 : msc.calculateLogP();
                logLikelihood = likelihood.getArrayValue();

                newLogLikelihood = logPriorProb + logMSC + (logLikelihood * beta); 

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
