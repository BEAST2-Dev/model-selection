package beast.gss;



import java.util.ArrayList;
import java.util.List;

import beast.app.util.Application;
import beast.core.Distribution;
import beast.core.*;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.cpo.BEASTRunAnalyser;
import beast.evolution.tree.TreeDistribution;
import beast.gss.distributions.CCProbability;
import beast.gss.distributions.KernelDensityEstimatorDistribution;
import beast.gss.distributions.MultivariateKDEDistribution;
import beast.gss.distributions.NormalKDEDistribution;
import beast.math.distributions.MRCAPrior;
import beast.util.LogAnalyser;

public class GeneralisedSteppingStone extends BEASTRunAnalyser {

	@Override
	public void initAndValidate() {
		//nothing to do
	}

	@Override
	public void run() throws Exception {
		MCMC mcmc = getMCMC();
		LogAnalyser tracelog = getTraceLog(mcmc);
		
		CompoundDistribution prior = getPrior(mcmc);
		CompoundDistribution altPrior = getAltPrior(prior, tracelog);
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
		return cd;
	}
	
	private Distribution getAltTreeDist(TreeDistribution d) {
		CCProbability ccDistr = new CCProbability(treeFileInput.get(), d.treeInput.get(), 10);
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
		int index = labels.indexOf(id);
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


	public static void main(String[] args) throws Exception {
		new Application(new GeneralisedSteppingStone(), "Generalised Stepping Stone", args);
	}
}
