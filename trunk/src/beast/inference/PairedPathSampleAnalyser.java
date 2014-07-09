package beast.inference;

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import beast.app.util.ConsoleApp;
import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.inference.PairedPathSampler.Scheme;
import beast.util.LogAnalyser;



@Description("Reads logs produces through PairedPathSampler and estimates Bayes factor of the pair of models")
public class PairedPathSampleAnalyser extends beast.core.Runnable {
	public Input<String> rootDirInput = new Input<String>("rootdir", "root directory for storing particle states and log files (default /tmp)", "/tmp");
	public Input<Double> alphaInput = new Input<Double>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
			"If alpha <= 0, uniform intervals are used.", 0.3);
	public Input<Integer> stepsInput = new Input<Integer>("nrOfSteps", "the number of steps to use, default 8", 8);
	public Input<Integer> burnInPercentageInput = new Input<Integer>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);

	DecimalFormat formatter;
	
	@Override
	public void initAndValidate() throws Exception {
	}
	
	/** estimate marginal likelihoods from logs produced by PathSampler
	 * @param nSteps number of steps used by PathSampler
	 * @param alpha  if < 0 uniform intervals are used, otherwise a Beta(alpha,1.0) distribution is used for intervals
	 * @param rootDir location where log files are stored
	 * @param burnInPercentage percentage of log files to be discarded
	 * @return log of marginal likelihood
	 * @throws Exception
	 */
	double estimateMarginalLikelihood(int nSteps, double alpha, String rootDir, int burnInPercentage) throws Exception {
		List<List<Double>> logdata = new ArrayList<List<Double>>(); 
		
		
		String sFormat = "";
		for (int i = nSteps; i > 0; i /= 10) {
			sFormat += "#";
		}
		formatter = new DecimalFormat(sFormat);

		// collect likelihood estimates for each step
		double [] marginalLs = new double[nSteps];
		Double [] [] marginalLs2 = new Double[nSteps][];
		for (int i = 0; i < nSteps; i++) {
			List<Double> logdata1 = new ArrayList<Double>();
			String logFile = getStepDir(rootDir, i) + "/" + PathSampler.LIKELIHOOD_LOG_FILE;
			LogAnalyser analyser = new LogAnalyser(new String[] {logFile}, 2000, burnInPercentage);
			marginalLs[i] = analyser.getMean("diff-posterior");
			marginalLs2[i] = analyser.getTrace("diff-posterior");
			System.err.println("marginalLs[" + i + " ] = " + marginalLs[i]);

			logdata1.add(marginalLs[i]);
			logdata1.add(0.0);
			logdata1.add(analyser.getESS("diff-posterior"));
			logdata.add(logdata1);
		}
		
		// combine steps
		double logBF = 0;
		PairedPathSampler.Scheme scheme = Scheme.uniform;
		if (alpha <= 0) { 
			scheme = Scheme.uniform;
			// uniform intervals
			for (int i = 0; i < nSteps - 1; i++) {
				logBF += (marginalLs[i] + marginalLs[i + 1]); 
			}		
			logBF = logBF / (2.0 * (nSteps - 1));
		} else {
			scheme = Scheme.sigmoid;
			// intervals follow Beta distribution
			//BetaDistribution betaDistribution = new BetaDistributionImpl(alpha, 1.0);
			double [] contrib = new double[nSteps-1];
			
			for (int i = 0; i < nSteps - 1; i++) {
				List<Double> logdata1 = logdata.get(i);;
				double beta1 = PairedPathSampler.nextBeta(scheme, i, nSteps - 1, alpha); 
				double beta2 = PairedPathSampler.nextBeta(scheme, i + 1,  nSteps - 1, alpha); 
				double weight = beta2 - beta1;

				// Use formula (18) 
				// Make the most of your samples: Bayes factor estimators for high-dimensional models of sequence evolution
				// G Baele, P Lemey, S Vansteelandt
				// BMC bioinformatics 14 (1), 85
				Double [] marginal2 = marginalLs2[i];
				double logLmax = max(marginal2);
				logBF += weight * logLmax;
				
				int n = marginal2.length;
				double x = 0;
				for (int j = 0; j < n; j++) {
					x += Math.exp(weight * (marginal2[j] - logLmax)); 
				}
				logBF += Math.log(x/n);

				contrib[i] = -(weight * logLmax + Math.log(x/n));
				logdata1.set(1, -(weight * logLmax + Math.log(x/n)));

//				logMarginalL += weight * marginalLs[i]; 
			}
						
		}
		
		System.out.println("\nStep         beta       " +
				"diff-posterior contribution ESS");
		//BetaDistribution betaDistribution = new BetaDistributionImpl(alpha, 1.0);
		for (int i = 0; i < nSteps; i++) {
			System.out.print(format(i)+" ");
			double beta = PairedPathSampler.nextBeta(scheme, i, nSteps-1, alpha); 
			System.out.print(format(beta)+" ");

			
			for (Double d : logdata.get(i)) {
				System.out.print(format(d) + " ");
			}
			System.out.println();
		}		
		System.out.println();
		return -logBF;
	}

	private String format(double d) {
		DecimalFormat format = new DecimalFormat("###.####");
		String s = format.format(d);
		if (s.length() < 12) {
			s += "            ".substring(s.length());
		}
		return s;
	}

	private double max(Double[] marginal2) {
		Double max = marginal2[0];
		for (Double v : marginal2) {
			max = Math.max(v, max);
		}
		return max;
	}

	String getStepDir(String rootDir, int iParticle) {
		return rootDir + "/step" + formatter.format(iParticle);
	}
	
	public static void main(String[] args) throws Exception {
		PairedPathSampleAnalyser analyser = new PairedPathSampleAnalyser();
		int nSteps = Integer.parseInt(args[0]);
		double alpha = Double.parseDouble(args[1]);
		String rootDir = args[2];
		int burnInPercentage = Integer.parseInt(args[3]);
		double marginalL = analyser.estimateMarginalLikelihood(nSteps, alpha, rootDir, burnInPercentage);
		System.out.println("Bayes factor estimate = " + marginalL);
	}
	
	
    ConsoleApp consoleApp = null;

	@Override
	public void run() throws IOException {
		// create output window
        String nameString = "PairedPathSampleAnalyser";
        String title = "Paired Path Sample Analyser -- " + rootDirInput.get();
        consoleApp = new ConsoleApp(nameString, title);

		double marginalL = Double.NaN;
		try {
			marginalL = estimateMarginalLikelihood(
					stepsInput.get(), 
					alphaInput.get(), 
					rootDirInput.get(), 
					burnInPercentageInput.get());
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Bayes factor estimate = " + marginalL);
	}
}
