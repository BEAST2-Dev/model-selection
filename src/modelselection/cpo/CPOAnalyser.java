package modelselection.cpo;



import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.treeannotator.TreeAnnotator;
import beastfx.app.treeannotator.TreeAnnotator.TreeSet;
import beastfx.app.methodsection.implementation.BranchRateModelMethodsText;
import beastfx.app.tools.Application;
import beastfx.app.util.LogFile;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.core.Param;
import beast.base.inference.Runnable;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.likelihood.ThreadedTreeLikelihood;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beastfx.app.tools.LogAnalyser;
import beast.base.util.Randomizer;
import beast.base.parser.XMLParser;
import beast.base.parser.XMLParserException;
import modelselection.core.CPOLogger;

@Description("Calculate Conditional Predictive Ordinates (CPO), which is a leave one out cross validation measure of fit "
		+ "as described in Lewis et al, Sys Bio, 2014, but adds bootstrap variance estimate as well.")	
// allows multiple partitions
// limitations: assumes a single tree (so does not work with *BEAST)
// not sure how to deal with continuous traits (e.g. geography)
public class CPOAnalyser extends BEASTRunAnalyser {

	final public Input<LogFile> cpoLogFileInput = new Input<>("cpologFile","input file containing CPO log produced by an earlier BEAST run using a CPOLogger -- only specify when XML is not specified");

	public Input<Integer> bootstrapLengthInput = new Input<>("bootstrapLength", "number of bootstrap samples used to calculate variance in CPO estimate" , 1000);

	double EPSILON = 1e-6; // may need an input for this
	
	
	public CPOAnalyser() {}
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		System.setProperty("java.only", "true");
		if (burninInput.get() < 0 || burninInput.get() > 100) {
			throw new IllegalArgumentException("burnin is a percentage and should not be less than 0 or larger than 100");
		}
		
		CPOTable cpoTable = null;
		if (isSpecified(xmlFileInput.get())) {
			cpoTable = getCPOTableFromXML();
			if (cpoTable == null) {
				cpoTable = runXMLToGetCPOTable();
			}
		} else {
			cpoTable = getCPOTableFromCPOLog(cpoLogFileInput.get());
		}
		
		double [][] patterLogProbs = cpoTable.patternLogProbs;
		int [] weights = cpoTable.patternWeights;
		int patternCount = patterLogProbs.length;
		int treeCount = patterLogProbs[0].length;

    	// precalc total log pseudomarginal likelihood (LPML) from siteProbs array
  		// Part of Equation (14) in Lewis et al 2014
    	double [] minLogP = new double[patternCount];
		for (int i = 0; i < patternCount; i++) {
			double min = patterLogProbs[i][0];
			for (int k = 0; k < treeCount; k++) {
    			min = Math.min(min, patterLogProbs[i][k]);
    		}
    		minLogP[i] = min;
    	}
    	
		int [] order = new int[treeCount];
		for (int i = 0; i < order.length; i++) {
			order[i] = i;
		}
		double LPML = calcLPML(order, minLogP, cpoTable);

  		Log.info("\nlog pseudomarginal likelihood (LPML) = " + LPML);
    	
    	Log.warning.println("Calculating variance of CPO");
    	Log.warning.println("|---------|---------|---------|---------|---------|---------|---------|---------|");
    	int siteCount = 0;
    	for (int d : weights) {
    		siteCount += d;
    	}
    	// maps site index to pattern index
    	int [] siteMap = new int[siteCount];
    	int k = 0;
    	for (int i = 0; i < patternCount; i++) {
    		for (int j = 0; j < weights[i]; j++) {
    			siteMap[k++] = i;
    		}
    	}
    	double reported = 0;
    	int replicates = bootstrapLengthInput.get();
    	Log.warning.println(replicates + " replicates");
    	k = 0;
    	//Randomizer.setSeed(1237);
    	double [] LPMLs = new double[replicates];
    	while (k < replicates) {
        	// calc CPO from subsample of siteProbs array
			order = Randomizer.sampleIndicesWithReplacement(treeCount);
    		LPMLs[k] = calcLPML(order, minLogP, cpoTable);
    		
			while (reported < k) {
				Log.warning.print("*");
				reported += replicates / 80;
			}
			k++;		
		}

    	// Summarise 
    	double mean = 0;
        for (int i = 0; i < replicates; i++) {
            mean += LPMLs[i];
        }
        mean /= replicates;
    	
        double var = 0;
        for (int i = 0; i < replicates; i++) {
            var += (LPMLs[i] - mean) *
                    (LPMLs[i] - mean);
        }
        var /= (replicates - 1.0);
        double standardDeviation = Math.sqrt(var);
        Log.info("\nmean = " + mean + ", standardDeviation = " + standardDeviation);

	}
		
	private CPOTable runXMLToGetCPOTable() throws IOException, SAXException, ParserConfigurationException, XMLParserException {
		XMLParser parser = new XMLParser();
		Runnable o = parser.parseFile(xmlFileInput.get());
		if (!(o instanceof MCMC)) {
			throw new IllegalArgumentException("XML should contain an MCMC analysis");
		}
		MCMC mcmc = (MCMC) o;
		mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
		File tmpFile = File.createTempFile("cpoLog","log");
		Logger traceLog = null;
		for (Logger logger : mcmc.loggersInput.get()) {
			if (logger.fileNameInput.get() != null && logger.modeInput.get() != Logger.LOGMODE.tree) {
				traceLog = logger;
				break;
			}
		}
		Logger cpoLogger = new CPOLogger();
		cpoLogger.initByName("fileName", tmpFile.getAbsolutePath(), 
				"logEvery", traceLog.everyInput.get(),
				"mode", traceLog.modeInput.get(),
				"sanitiseHeaders", traceLog.sanitiseHeadersInput
				);
		mcmc.loggersInput.setValue(cpoLogger, mcmc);
		mcmc.initAndValidate();
		mcmc.run();
		
		CPOTable table = getCPOTableFromCPOLog(tmpFile);
		return table;
	}

	private CPOTable getCPOTableFromCPOLog(File cpoLogFile) throws IOException {
		CPOTable table = new CPOTable();
		table.readFromCPOLog(cpoLogFile);
		return table;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public CPOTable getCPOTableFromXML() throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		// collect data
		MCMC mcmc = getMCMC();

		LogAnalyser tracelog = getTraceLog(mcmc, traceFileInput.get(), burninInput.get());

		TreeAnnotator.TreeSet treeSet = getTreeSet(mcmc);
    	treeSet.reset();
    	int treeSetSize = treeSetSize(treeSet);

		// sanity check
    	if (treeSetSize != tracelog.getTrace(0).length) {
    		throw new IllegalArgumentException("Trace and log files appear to be sampled with different frequencies (cannot hanlde that yet)");
    	}
    	
    	Double [] likelihoods = tracelog.getTrace("likelihood");
    	
    	// reserve memory for site probabilities
    	CompoundDistribution likelihood = getLikelihood(mcmc); 
    	int patternCount = 0;
    	for (Distribution d : likelihood.pDistributions.get()) {
    		if (d instanceof GenericTreeLikelihood) {
    			patternCount += ((GenericTreeLikelihood) d).dataInput.get().getPatternCount();
    		}
    	}
    	
    	CPOTable cpoTable = new CPOTable();
    	cpoTable.patternLogProbs = new double[patternCount][treeSetSize];
    	cpoTable.patternWeights = new int[patternCount];
    	double [][] patterLogProbs = cpoTable.patternLogProbs;
    	int k = 0;
    	for (Distribution d : likelihood.pDistributions.get()) {
    		if (d instanceof GenericTreeLikelihood) {
    			int [] currentWeights = ((GenericTreeLikelihood) d).dataInput.get().getWeights();
    			System.arraycopy(currentWeights, 0, cpoTable.patternWeights, k, currentWeights.length);
    			k += currentWeights.length;
    		}
    	}

    	// process State
    	State state = mcmc.startStateInput.get();
    	List<StateNode> stateNodes = state.stateNodeInput.get();
    	List<Object> values = new ArrayList<>();
    	List<String> labels = tracelog.getLabels();
    	Tree tree = null;
    	for (int i = 0; i < stateNodes.size(); i++) {
    		StateNode stateNode = stateNodes.get(i);
    		String id = stateNode.getID();
    		String shortid = id.contains(".") ? id.substring(0, id.lastIndexOf('.')): id;
    		if (stateNode instanceof Parameter<?> && ((Parameter<?>) stateNode).getDimension() > 1) {
    			Parameter<?> p = ((Parameter<?>) stateNode);
    			id += ".1";
    			int index = labels.indexOf(shortid);
    			if (index < 0) {
        			shortid += ".1";
    				index = labels.indexOf(shortid);
    			}
    			if (index >= 0) {
    				Double [][] _values = new Double[p.getDimension()][];
    				for (int j = 0; j < p.getDimension(); j++) {
    					Double [] v = tracelog.getTrace(1 + index + j);
    					_values[j] = v;
    				}
    				values.add(_values);
    			} else {
    				Log.warning.println("Did not find " + id + " in log.");
    				values.add(null);
    			}
    		} else if (stateNodes.get(i) instanceof Tree) {
    			tree = (Tree) stateNodes.get(i);
    			values.add(null);
    		} else {
    			int index = labels.indexOf(id);
    			if (index < 0) {
    				index = labels.indexOf(shortid);
    			}
    			if (index >= 0) {
    				Double [] v = tracelog.getTrace(index + 1);
    				values.add(v);
    			} else {
    				Log.warning.println("Did not find " + id + " in log.");
    				values.add(null);
    			}
    		}
    	}
    	
    	setUpBranchRateModel(likelihood);
    	    	
		Log.warning.println("Calculating total CPO");
    	Log.warning.println("|---------|---------|---------|---------|---------|---------|---------|---------|");
    	k = 0;
		int reported = 0;
		treeSet.reset();
    	while (treeSet.hasNext()) {
			// set up the state
    		for (int i = 0; i < stateNodes.size(); i++) {
    			if (values.get(i) != null) {
    				Object o = values.get(i);
					StateNode stateNode = stateNodes.get(i);
    				if (o instanceof Double[]) {
    					double value = ((Double[])o)[k];
						((Parameter)stateNode).setValue(value);
    				} else {
    					Double [][] _values = (Double[][]) o;
    					for (int j = 0; j < _values.length; j++) {
    						((Parameter)stateNode).setValue(j, _values[j][k]);
    					}
    				}
    			}
    		}
			Tree currentTree = treeSet.next();
			//TreeParser p = new TreeParser(tree.getTaxonset().asStringList(), newick, 0, false);
			tree.assignFrom(currentTree);

			// calc likelihood
			mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
			// grab site probabilities from treelikelihoods
			int i = 0;
	    	for (Distribution d : likelihood.pDistributions.get()) {
	    		if (d instanceof TreeLikelihood) {
	    			TreeLikelihood tl = (TreeLikelihood) d;
	    			double [] pll = tl.getPatternLogLikelihoods();
	    			for (int j = 0; j < pll.length; j++) {
	    				patterLogProbs[i++][k] = pll[j];
	    			}
	    		} else if (d instanceof ThreadedTreeLikelihood) {
	    			ThreadedTreeLikelihood tl = (ThreadedTreeLikelihood) d;
	    			double [] pll = tl.getPatternLogLikelihoods();
	    			for (int j = 0; j < pll.length; j++) {
	    				patterLogProbs[i++][k] = pll[j];
	    			}
	    		}

	    	}
			double lastLikelihood  = likelihood.getCurrentLogP();
			if (Math.abs(lastLikelihood - likelihoods[k]) > EPSILON) {
				Log.warning("The difference between calculated likelihood and likelihood in log file is " + (lastLikelihood - likelihoods[k]) + 
						" larger than " + EPSILON + ". This indicates that the state cannot be restored reliably from the log files." +
						" Giving up site probablity reconstruction from XML & log files.");
				return null;
			}
					
			while (reported - k < 0) {
				Log.warning.print("*");
				reported += 1 + treeSetSize / 86;
			}
			k++;		
		}

    	return cpoTable;
	}



	private int treeSetSize(TreeSet treeSet) throws IOException {
		int size = 0;
		treeSet.reset();
    	while (treeSet.hasNext()) {
    		treeSet.next();
    		size++;
    	}
		treeSet.reset();
		return size;
	}

	private void setUpBranchRateModel(CompoundDistribution likelihood) {
    	// set up branch rate model to pick up rates from Tree metadata
		for (Distribution d : likelihood.pDistributions.get()) {
			if (d instanceof GenericTreeLikelihood) {
				GenericTreeLikelihood tl = (GenericTreeLikelihood) d;
				BEASTInterface clock = (BEASTInterface) tl.branchRateModelInput.get();
				if (hasInput(clock, "tree")) {
					Input<?> input = clock.getInput("tree");
					if (input != null) {
						Object o = input.get();
						if (o != null &&  o instanceof Tree) {
							Tree tree = (Tree) o;
							RateByMetaData newClock = new RateByMetaData(tree);
							newClock.meanRateInput.setValue(clock.getInput("clock.rate").get(), newClock);
							tl.branchRateModelInput.setValue(newClock, tl);
							tl.initAndValidate();
						}
					}
				}
			}
		}
	}

	private boolean hasInput(BEASTInterface clock, String name) {
		for (Input<?> input : clock.listInputs()) {
			if (input.getName().equals(name)) {
				return true;
			}
		}
		return false;
	}

	@Description("clock model that gets its rates from metadata on a newick tree")
	class RateByMetaData extends BranchRateModel.Base {

		RateByMetaData(@Param(name="tree",description="beast tree with metadata containing rates")Tree tree) {
			this.tree = tree;
		}
		
		@Override
		public void initAndValidate() {
		}

		@Override
		public double getRateForBranch(Node node) {
			Node src = tree.getNode(node.getNr());
			Object o = src.getMetaData("rate");
			if (o == null) {
				return meanRateInput.get().getArrayValue();
			}
			if (o instanceof Double) {
				return (Double) o;
			}
			String value = o.toString();
			double rate = Double.parseDouble(value);
			return rate;
		}

		Tree tree;
		public Tree getTree() {return tree;}
		public void setTree(Tree tree) {this.tree = tree;}
	}
	
	@Description("Class containing CPO read from a log file")
	public
	class CPOTable {
		double [][] patternLogProbs;
		int [] patternWeights;
		
		
		// expects file to be a file where the first line contain patternweights in this format:
		// #weights 32 12 32 1 3 4 5 ...
		// and the remainder a standard BEAST log file
		void readFromCPOLog(File file) throws IOException {
			LogAnalyser cpoLog = new LogAnalyser(file.getAbsolutePath(), burninInput.get(), false, false);
			patternLogProbs = new double [cpoLog.getLabels().size()][cpoLog.getTrace(0).length];
			for (int i = 0; i < patternLogProbs.length; i++) {
				Double [] d = cpoLog.getTrace(i + 1);
				for (int j = 0; j < d.length; j++) {
					patternLogProbs[i][j] = d[j];
				}
			}
			
			patternWeights = new int[patternLogProbs.length];
	        BufferedReader fin = new BufferedReader(new FileReader(file));
	        String str = fin.readLine();
        	fin.close();
	        String [] strs = str.split("\\s+");
	        if (strs.length != patternWeights.length + 1) {
	        	throw new IllegalArgumentException("weights in file are not equal to columns in file");
	        }
	        for (int i = 0; i < patternWeights.length; i++) {
	        	patternWeights[i] = Integer.parseInt(strs[i+1]);
	        }
	    }
		
		@Override
		public boolean equals(Object o) {
			if (!(o instanceof CPOTable)) {
				return false;
			}
			CPOTable table2 = (CPOTable) o;
			if (patternWeights.length != table2.patternWeights.length) {
				return false;
			}
			for (int i = 0; i < patternWeights.length; i++) {
				if (patternWeights[i] != table2.patternWeights[i]) {
					return false;
				}
			}
			if (patternLogProbs.length != table2.patternLogProbs.length ||
				patternLogProbs[0].length != table2.patternLogProbs[0].length) {
				return false;
			}
			for (int i = 0; i < patternLogProbs.length; i++) {
				for (int j = 0; j < patternLogProbs[0].length; j++) {
					if (Math.abs(patternLogProbs[i][j] - table2.patternLogProbs[i][j]) > 1e-10) {
						return false;
					}
 				}
			}
			return true;
		}

		public double[][] getPatternLogProbs() {
			return patternLogProbs;
		}

		public int[] getPatternWeights() {
			return patternWeights;
		}
	}
	
	
	
	public double calcLPML(int [] order, double[] minLogP, CPOTable cpoTable) {
		double [][] patterLogProbs = cpoTable.patternLogProbs;
		int patternCount = patterLogProbs.length;
		int treeSetSize = patterLogProbs[0].length;
		int [] patternWeights = cpoTable.patternWeights;
  		// Equation (14) in Lewis et al 2014
    	double [] CPOi = new double[patternCount];
  		for (int i = 0; i < patternCount; i++) {
    		double CPOi_ = Math.log(treeSetSize) + minLogP[i];
    		double sum = 0;
    		double [] p = patterLogProbs[i];
        	for (int k : order) {
        		sum += Math.exp(minLogP[i] - p[k]);
    		}
        	CPOi_ -= Math.log(sum);
    		CPOi[i] = CPOi_;
    	}
    	
  		// Equation (15) in Lewis et al 2014
  		double LPML = 0;
  		for (int i = 0; i < patternCount; i++) {
  			LPML += CPOi[i] * patternWeights[i];
  		}
		return LPML;
	}

	private CompoundDistribution getLikelihood(MCMC mcmc) {
		if (!(mcmc.posteriorInput.get() instanceof CompoundDistribution)) {
			throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with posterior CompoundDistribution");
		}
		CompoundDistribution posterior = (CompoundDistribution) mcmc.posteriorInput.get();
		for (Distribution d : posterior.pDistributions.get()) {
			if (d.getID().equals("likelihood")) {
				if (!(d instanceof CompoundDistribution)) {
					throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with likelihood CompoundDistribution");
				}
				return (CompoundDistribution) d;
			}
		}
		throw new IllegalArgumentException("The XML file does not seem to contain an MCMC analysis with a 'likelihood' in the posterior");
	} // getLikelihood

	public static void main(String[] args) throws Exception {
		new Application(new CPOAnalyser(), "CPOAnalyser", args);
	}
}



