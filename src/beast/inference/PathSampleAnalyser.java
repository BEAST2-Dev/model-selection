package beast.inference;


import beast.app.util.Application;
import beast.app.util.ConsoleApp;
import beast.core.Description;
import beast.core.Input;
import beast.util.HeapSort;
import beast.util.LogAnalyser;
import beast.util.Randomizer;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.BetaDistribution;
import org.apache.commons.math.distribution.BetaDistributionImpl;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;



@Description("Reads logs produces through PathSampler and estimates marginal likelihood")
public class PathSampleAnalyser extends beast.core.Runnable {
	public Input<String> rootDirInput = new Input<>("rootdir", "root directory for storing particle states and log files (default /tmp)", "/tmp");
	public Input<Double> alphaInput = new Input<>("alpha", "alpha parameter of Beta(alpha,1) distribution used to space out steps, default 0.3" +
			"If alpha <= 0, uniform intervals are used.", 0.3);
	public Input<Integer> stepsInput = new Input<>("nrOfSteps", "the number of steps to use, default 8", 8);
	public Input<Integer> burnInPercentageInput = new Input<>("burnInPercentage", "burn-In Percentage used for analysing log files", 50);
	public Input<Integer> crossValInput = new Input<>("cross", "the number of cross validation intervals to use for estimating the variance, default 0, no variance estimate is provided", 0);
	public Input<Integer> crossVaRepeatslInput = new Input<>("repeats", "the number of times cross validation is repeated, default 100", 100);
	public Input<Boolean> getBetaFromFileInput = new Input<>("betaFromFile", "get beta values from XML file instead of assuming they are from a Beta distirbution", false);
	public Input<String> tagInput = new Input<>("tag", "label of the column in the likelihood.log file", "likelihood");
	
	DecimalFormat formatter;
	
	@Override
	public void initAndValidate() {
	}
	
	/** estimate marginal likelihoods from logs produced by PathSampler
	 * @param nSteps number of steps used by PathSampler
	 * @param alpha  if < 0 uniform intervals are used, otherwise a Beta(alpha,1.0) distribution is used for intervals
	 * @param rootDir location where log files are stored
	 * @param burnInPercentage percentage of log files to be discarded
	 * @return log of marginal likelihood
	 * @throws Exception
	 */
	public double estimateMarginalLikelihood(int nSteps, double alpha, String rootDir, int burnInPercentage) throws Exception {
		return estimateMarginalLikelihood(nSteps, alpha, rootDir, burnInPercentage, "likelihood");
	}
	
	public double estimateMarginalLikelihood(int nSteps, double alpha, String rootDir, int burnInPercentage, String tag) throws Exception {
		String sFormat = "";
		for (int i = nSteps; i > 0; i /= 10) {
			sFormat += "#";
		}
		formatter = new DecimalFormat(sFormat);

		List<List<Double>> logdata = new ArrayList<List<Double>>(); 
		
		
		double [] beta = new double[nSteps];
		if (getBetaFromFileInput.get()) {
			for (int i = 0; i < nSteps; i++) {
				beta[i] = 1.0;
				String xmlFile = getStepDir(rootDir, i) + "/beast.xml";
		        BufferedReader fin = new BufferedReader(new FileReader(xmlFile));
		        while (fin.ready()) {
		            String str = fin.readLine();
		            if (str.matches(".*PathSamplingStep.*")) {
		            	if (str.matches(".*beta=['\"].*")) {
		            		str = str.replaceAll(".*beta=['\"]", "");
		            		str = str.replaceAll("['\"].*", "");
		    				beta[i] = Double.parseDouble(str);
		            	}
		            }
		        }
		        fin.close();
			}
		} else {
			BetaDistribution betaDistribution = null;
			if (alphaInput.get() > 0){
				betaDistribution = new BetaDistributionImpl(alphaInput.get(), 1.0);
			}

			for (int i = 0; i < nSteps; i++) {
				beta[i] = betaDistribution != null ?
						betaDistribution.inverseCumulativeProbability((nSteps - 1.0 - i)/ (nSteps - 1))
						:(nSteps - 1.0 - i)/ (nSteps - 1);

			}
		}
		
		int [] order = new int[nSteps];
        HeapSort.sort(beta, order);
        Arrays.sort(beta);

		// collect likelihood estimates for each step
		double [] marginalLs = new double[nSteps];
		Double [] [] marginalLs2 = new Double[nSteps][];
		for (int i = 0; i < nSteps; i++) {
			List<Double> logdata1 = new ArrayList<Double>();
			String logFile = getStepDir(rootDir, i) + "/" + PathSampler.LIKELIHOOD_LOG_FILE;
			LogAnalyser analyser = new LogAnalyser(new String[] {logFile}, burnInPercentage);
			marginalLs[order[i]] = analyser.getMean(tag);
			marginalLs2[order[i]] = analyser.getTrace(tag);
			System.out.println("marginalLs[" + i + " ] = " + marginalLs[order[i]]);

			logdata1.add(marginalLs[order[i]]);
			logdata1.add(0.0);
			logdata1.add(analyser.getESS(tag));
			logdata.add(logdata1);
		}
		double logMarginalL = estimateMarginalLikelihood(logdata, marginalLs2, marginalLs, alpha, nSteps, beta, true);

		Double [] [] subMarginalLs = new Double[nSteps][];
		int n = crossValInput.get();
		if (n > 1) {
			// calculate ML for each of the cross validated subsets 
			for (int i = 0; i < nSteps; i++) {
				subMarginalLs[i] = new Double[marginalLs2[i].length * (n - 1) / n];
			}
			int REPEATS = crossVaRepeatslInput.get();
			
			double [][] m = new double[n][REPEATS];
			
			for (int repeats = 0; repeats < REPEATS; repeats++) {
				for (int j = 0; j < nSteps; j++) {
					randomise(marginalLs2[j]);
				}				
				
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < nSteps; j++) {
						int lo = i * marginalLs2[j].length / n;
						int hi = i+1 * marginalLs2[j].length / n;
						int k = 0;
						int r = 0;
						while (k < lo) {
							subMarginalLs[j][r] = marginalLs2[j][k];
							k++;
							r++;
						}
						k = hi;
						while (r < subMarginalLs[j].length) {
							subMarginalLs[j][r] = marginalLs2[j][k];
							k++;
							r++;
						}
						marginalLs[j] = mean(subMarginalLs[j]);
					}
					m[i][repeats] = estimateMarginalLikelihood(logdata, subMarginalLs, marginalLs, alpha, nSteps, beta, false);
				}
				System.err.print('.');
			}
			System.err.println();
			
			double mean = 0;
			for (int repeats = 0; repeats < REPEATS; repeats++) {
				for (int i = 0; i < n; i++) {
					mean += m[i][repeats];
				}
			}
			mean /= (n * REPEATS);
			
			double var = 0;
			for (int repeats = 0; repeats < REPEATS; repeats++) {
				for (int i = 0; i < n; i++) {
					var += (m[i][repeats] - mean ) * (m[i][repeats] - mean);
				}
			}
			var /= (n * REPEATS - 1);
			double sd = Math.sqrt(var);
			// factor compensating for correlation between samples
			double FUDGE = 2.0;
			System.out.println("SD: " + sd * FUDGE);
		}
		
		return logMarginalL;
	}
	
	private void randomise(Double[] doubles) {
		int n = doubles.length;
		for (int i = 0; i < n; i++) {
			int t = Randomizer.nextInt(n);
			double tmp = doubles[i];
			doubles[i] = doubles[t];
			doubles[t] = tmp;
		}
		
	}

	private double mean(Double[] doubles) {
		double sum = 0;
		for (Double d : doubles) {
			if (d == null) {
				int h = 3;
				h++;
			}
			sum += d;
		}
		sum /= doubles.length;
		return sum;
	}

	private double estimateMarginalLikelihood(List<List<Double>> logdata, Double [] [] marginalLs2, double [] marginalLs, double alpha, int nSteps, double [] betas, boolean verbose) throws MathException, InterruptedException {
		// combine steps
		double logMarginalL = 0;
		if (alpha <= 0) { 
			// uniform intervals
			for (int i = 0; i < nSteps - 1; i++) {
				logMarginalL += (marginalLs[i] + marginalLs[i + 1]); 
				logdata.get(i).set(1, (marginalLs[i] + marginalLs[i + 1])/(2.0 * (nSteps - 1)));
			}		
			logMarginalL = logMarginalL / (2.0 * (nSteps - 1));
		} else {
			// intervals follow Beta distribution
			double [] contrib = new double[nSteps-1];
			
			for (int i = 0; i < nSteps - 1; i++) {
				List<Double> logdata1 = logdata.get(nSteps-i-1);
				double beta1 = betas[i]; 
				double beta2 = betas[i+1]; 
				double weight = beta2 - beta1;

				// Use formula top right at page 153 of 
				// Xie W, Lewis PO, Fan Y, Kuo L, Chen MH. 2011. Improving marginal
				// likelihood estimation for Bayesian phylogenetic model selection.
				// Syst Biol. 60:150-160.
				Double [] marginal2 = marginalLs2[i];
				double logLmax = max(marginal2);
				logMarginalL += weight * logLmax;
				
				int n = marginal2.length;
				double x = 0;
				for (int j = 0; j < n; j++) {
					x += Math.exp(weight * (marginal2[j] - logLmax)); 
				}
				logMarginalL += Math.log(x/n);

				contrib[i] = (weight * logLmax + Math.log(x/n));
				logdata1.set(1, (weight * logLmax + Math.log(x/n)));

//				logMarginalL += weight * marginalLs[i]; 
			}
						
		}
		
		if (verbose) {
			if (consoleApp != null) {
				// allow output to flush to app window
				Thread.sleep(500);
			}
	
			System.out.println("\nStep        theta         likelihood   contribution ESS");
			BetaDistribution betaDistribution = alpha > 0 ? new BetaDistributionImpl(alpha, 1.0) : null;
			for (int i = 0; i < nSteps; i++) {
				System.out.print(format(i)+" ");
				double beta = betaDistribution != null ?
						betaDistribution.inverseCumulativeProbability((nSteps - 1.0 - i)/ (nSteps - 1)):
							(nSteps - 1.0 - i)/ (nSteps - 1);
				System.out.print(format(beta)+" ");
	
				
				for (Double d : logdata.get(i)) {
					System.out.print(format(d) + " ");
				}
				System.out.println();
			}		
			double sumESS = 0;
			for (int i = 0; i < nSteps; i++) {
				sumESS += logdata.get(i).get(2);
			}
			System.out.println("sum(ESS) = " + format(sumESS));
			System.out.println();
		}
		return logMarginalL;
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
		new Application(new PathSampleAnalyser(), "Path Sample Analyser", args);
		
//		PathSampleAnalyser analyser = new PathSampleAnalyser();
//		int nSteps = Integer.parseInt(args[0]);
//		double alpha = Double.parseDouble(args[1]);
//		String rootDir = args[2];
//		int burnInPercentage = Integer.parseInt(args[3]);
//		double marginalL = analyser.estimateMarginalLikelihood(nSteps, alpha, rootDir, burnInPercentage);
//		System.out.println("marginal L estimate = " + marginalL);
	}

    public ConsoleApp consoleApp = null;
    
	@Override
	public void run() throws Exception {
		// create output window
//        String nameString = "PathSampleAnalyser";
//        String title = "Path Sample Analyser -- " + rootDirInput.get();
//        consoleApp = new ConsoleApp(nameString, title);
//        
        // do the work
        double marginalL = estimateMarginalLikelihood(
				stepsInput.get(), 
				alphaInput.get(), 
				rootDirInput.get(), 
				burnInPercentageInput.get(),
				tagInput.get());
        
        
		Thread.sleep(500);
		System.out.println("marginal L estimate = " + marginalL);
	}
	
}
