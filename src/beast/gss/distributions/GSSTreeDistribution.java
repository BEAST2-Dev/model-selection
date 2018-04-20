package beast.gss.distributions;


import java.io.IOException;
import java.util.*;

import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.TreeSet;
import beast.app.util.TreeFile;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.coalescent.TreeIntervals;

@Description("Tree Distribution consisting of a "
		+ "Conditional Clade Distribution for a tree set (as defined in "
		+ "Larget B (2013) The estimation of tree posterior probabilities "
		+ "using conditional clade probability distributions. Systematic Biology "
		+ "62(4): 501-511. http://dx.doi.org/10.1093/sysbio/syt014) and "
		+ "tree interval distribution.")
public class GSSTreeDistribution extends Distribution {
	public enum BranchLengthDistribution {useGamma, useIntervals, none};
	
	public Input<TreeFile> treeFileInput = new Input<>("treefile", "file containing tree set");
	public Input<TreeInterface> treeInput = new Input<>("tree", "beast tree for which the conditional clade distribution is calculated");
	public Input<Integer> burninPercentageInput = new Input<>("burnin", "percentage of the tree set to remove from the beginning", 10);
	public Input<BranchLengthDistribution> useGammaForBranchLengthsInput = new Input<>("useGammaForBranchLengths", "use an empirical gamma distribution for branch length distribution", BranchLengthDistribution.none, BranchLengthDistribution.values());

	
	

	private TreeFile treeFile;
	private TreeInterface tree;
	private Integer burninPercentage;
	private BranchLengthDistribution useGammaForBranchLengths = BranchLengthDistribution.none;
	
	private Tree lastTree;
	public Tree getLastTre() {return lastTree;}
	
	public TreeFile getTreefile() {return treeFile;}
	public void setTreefile(TreeFile treeFile) {this.treeFile = treeFile;}
	public TreeInterface getTree() {return tree;}
	public void setTree(TreeInterface tree) {this.tree = tree;}
	
    public Integer getBurnin() {
		return burninPercentage;
	}
	public void setBurnin(Integer burninPercentage) {
		this.burninPercentage = burninPercentage;
	}

	protected Map<BitSet, Map<BitSet, Clade>> conditionalCladeMap = new HashMap<>();
    protected Map<BitSet, Clade> cladeMap = new LinkedHashMap<>();
    
    private Map<String, Integer> mapTaxonIDToNr = new LinkedHashMap<>();
    
    NormalKDEDistribution [] distrs; 
    GammaDistribution gammaDistr;

    // log probability for a clade that does not exist in the clade system
    final static double EPSILON = -1e8;
	
    public GSSTreeDistribution() {}
//	public GSSTreeDistribution(@Param(name="treefile", description="file containing tree set") TreeFile treeFile,
//			@Param(name="tree", description="beast tree for which the conditional clade distribution is calculated") TreeInterface tree,
//			@Param(name="burnin", description="percentage of the tree set to remove from the beginning") Integer burninPercentage,
//			@Param(name="useGammaForBranchLengths", description="use an empirical gamma distribution for branch length distribution", defaultValue="none", optional=true) BranchLengthDistribution useGammaForBranchLengths) {
//		this.treeFile = treeFile;
//		this.tree = tree;
//		this.burninPercentage = burninPercentage;
//		if (burninPercentage < 0 || burninPercentage >= 100) {
//			throw new IllegalArgumentException("burnin must be a positive number not larger than 100");
//		}
//		this.useGammaForBranchLengths = useGammaForBranchLengths;
//		processTreeFile();
//	}
	
	
	@Override
	public void initAndValidate() {
		this.treeFile = treeFileInput.get();
		this.tree = treeInput.get();
		
		for (Node n : tree.getNodesAsArray()) {
			if (n.isLeaf()) {
				mapTaxonIDToNr.put(n.getID(), n.getNr());
			}
		}
				
		this.burninPercentage = burninPercentageInput.get();
		if (burninPercentage < 0 || burninPercentage >= 100) {
			throw new IllegalArgumentException("burnin must be a positive number not larger than 100");
		}
		this.useGammaForBranchLengths = useGammaForBranchLengthsInput.get();
		processTreeFile();
	}
	
	private void processTreeFile() {		
		try {

			TreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.getAbsolutePath(), burninPercentage);
			trees.reset();

			List<List<Double>> intervalLog = new ArrayList<>();
			Tree tree = trees.next();
			trees.reset();
    		for (int i = 0; i < tree.getNodeCount() + 1; i++) {
    			intervalLog.add(new ArrayList<>());
    		}
    		lastTree = tree;

			while (trees.hasNext()) {
				tree = trees.next();
				addClades(tree.getRoot());
				
				switch (useGammaForBranchLengths) {
				case useGamma:
					updateGammaEstimates(tree);
					break;
				case useIntervals:
					addToIntervalLog(tree, intervalLog);
					intervalLog.get(intervalLog.size() - 1).add(tree.getRoot().getHeight());
					break;
				case none:
				}
			}

//			for(BitSet b : cladeMap.keySet()) {
//				System.out.println(b + " " + cladeMap.get(b).count);
//			}
			
			if (this.tree == null) {
				this.tree = tree;
			}			
			
			switch (useGammaForBranchLengths) {
			case useGamma:
				gammaDistr = createGammaDistr();
				break;
			case useIntervals:
				createIntervalDistr(intervalLog);
				break;
			case none:
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
	private GammaDistribution createGammaDistr() {
		meanBranchLen = meanBranchLen / branchLenCount;
		double s = Math.log(meanBranchLen) - branchLogLen / branchLenCount;
		double alpha = 3 - s + Math.sqrt((s-3)*(s-3)+24*s)/(12*s); 
		double beta = meanBranchLen / alpha;
		
		GammaDistribution distr = new GammaDistributionImpl(alpha, beta);
		return distr;
	}

	public BranchLengthDistribution getUseGammaForBranchLengths() {return useGammaForBranchLengths;}
	public void setUseGammaForBranchLengths(BranchLengthDistribution useGammaForBranchLengths) {this.useGammaForBranchLengths = useGammaForBranchLengths;}
	public int getBranchLenCount() {return branchLenCount;}
	public void setBranchLenCount(int branchLenCount) {this.branchLenCount = branchLenCount;}

	double meanBranchLen = 0;
	int branchLenCount = 0;
	double branchLogLen = 0;
	
    private void updateGammaEstimates(Tree tree) {
    	for (Node node : tree.getNodesAsArray()) {
    		if (!node.isRoot()) {
    			double len = node.getLength();
    			meanBranchLen += len;
    			branchLogLen += Math.log(len);
    			branchLenCount++;
    		}
    	}
		
	}
	private void createIntervalDistr(List<List<Double>> intervalLog) {
		NormalKDEDistribution [] distrs = new NormalKDEDistribution[intervalLog.size()];
		for (int i = 0; i < distrs.length; i++) {
			List<Double> current = intervalLog.get(i);
			int n = current.size();
			Double [] sample = new Double[n];
			for (int j = 0; j < n; j++) {
				sample[j] = current.get(j); 
			}
			distrs[i] = new NormalKDEDistribution(sample, null);
		}
		this.distrs = distrs;
	}
    
	private void addToIntervalLog(Tree tree, List<List<Double>> intervalLog) {
    	TreeIntervals intervals = new TreeIntervals(tree);
    	for (int i = 0; i < intervals.getIntervalCount(); i++) {
    		List<Double> heights = intervalLog.get(i);
    		heights.add(intervals.getIntervalTime(i));
    	}
	}

    private BitSet addClades(Node node) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = mapTaxonIDToNr.get(node.getID());
            bits.set(index);

        } else {

        	BitSet left = addClades(node.getLeft());
        	BitSet right = addClades(node.getRight());
            bits.or(left);
            bits.or(right);

//            for (int i=1; i<bits.length(); i=i+2) {
//                bits.set(i, false);
//            }
            addClade(bits, left, right);
        }

        return bits;
    }

    private void addClade(BitSet bits, BitSet left, BitSet right) {
    	Clade clade = cladeMap.get(bits);

        if (clade == null) {
            clade = new Clade(bits);
            cladeMap.put(bits, clade);
        }
        clade.setCount(clade.getCount() + 1);
        	
        Map<BitSet, Clade> clades = conditionalCladeMap.get(bits);
        if (clades == null) {
        	clades = new HashMap<>();
            conditionalCladeMap.put(bits, clades);
        }
        
        clade = clades.get(left);
        if (clade == null) {
        	clade = new Clade(left);
        	clades.put(left, clade);
        }
        clade.setCount(clade.getCount() + 1);
        
        clade = clades.get(right);
        if (clade == null) {
        	clade = new Clade(right);
        	clades.put(right, clade);
        }
        clade.setCount(clade.getCount() + 1);
    }
	
	@Override
	public double calculateLogP() {
		logP = getLogCladeCredibility(tree.getRoot(), new BitSet());
		//System.err.print("GSST logp = " + logP);
		switch (useGammaForBranchLengths) {
		case useGamma:
			logP += getGammaBranchLengths(tree);
			break;
		case useIntervals:
			logP += getLogIntervalProbability(tree);
			break;
		case none:
			break;
		}
		//System.err.println(" GSST logp = " + logP);
		return logP;
	}
	
    private double getGammaBranchLengths(TreeInterface tree) {
    	 double logP = 0;
    	 for (Node node : tree.getNodesAsArray()) {
    		 double len = node.getLength();
    		 if (len > 0) {
    			 logP += gammaDistr.logDensity(len);
    		 }
    	 }
		return logP;
	}
	private double getLogIntervalProbability(TreeInterface tree) {
		TreeIntervals intervals = new TreeIntervals((Tree) tree);
		double logP = 0;
		for (int i = 0; i < intervals.getIntervalCount(); i++) {
			double t = intervals.getIntervalTime(i);
			if (t > 0) {
				double logPdf = distrs[i].logPdf(t);
				if (Double.isInfinite(logPdf)) {
					logPdf = EPSILON;
				}
				logP += logPdf;
			}
		}
		double logPdf = distrs[distrs.length - 1].logPdf(tree.getRoot().getHeight());
		if (Double.isInfinite(logPdf)) {
			logPdf = EPSILON;
		}
		logP += logPdf;
		return logP;
	}
	
//	class BitSet {
//		boolean [] b;
//		BitSet() {
//			b = new boolean[tree.getNodeCount()];
//		}
//		boolean get(int i) {
//			return b[i];
//		}
//		void set(int i) {
//			b[i] = true;
//		}
//		void set(int i, boolean value) {
//			b[i] = value;
//		}
//		void or(BitSet other) {
//			for (int i = 0;i < b.length; i++) {
//				b[i] |= other.b[i];
//			}
//		}
//		void xor(BitSet other) {
//			for (int i = 0;i < b.length; i++) {
//				if (!b[i] && other.b[i]) {
//					b[i] = true;
//				} else if (b[i] && !other.b[i]) {
//					b[i] = true;
//				} else {
//					b[i] = false;
//				}
//			}
//		}
//		int length() {
//			return b.length;
//		}
//		
//		@Override
//		public String toString() {
//			StringBuilder buf = new StringBuilder();
//			buf.append("{");
//			for (int i = 0; i < b.length; i++) {
//				if (b[i]) {
//					buf.append(i);
//					buf.append(' ');
//					buf.append(tree.getNode(i).getID());
//					buf.append(", ");
//				}
//			}
//			if (buf.length() > 1) {
//				buf.delete(buf.length()-2, buf.length()-1);
//			}
//			buf.append("}");					
//			return buf.toString();
//		}
//		
//		@Override
//		public boolean equals(Object obj) {
//			if (obj == null || !(obj instanceof BitSet)) {
//				return false;
//			}
//			BitSet other = (BitSet) obj;
//			if (other.b.length != b.length) {
//				return false;
//			}
//			if (obj == this) {
//				return true;
//			}
//			for (int i = 0; i < b.length; i++) {
//				if (b[i] != other.b[i]) {
//					return false;
//				}
//			}
// 			return true;
//		}
//	}
	
	public double getLogCladeCredibility(Node node, BitSet bits) {

        double logCladeCredibility = 0.0;

        if (node.isLeaf()) {
            int index = node.getNr();
            bits.set(index);
        } else {

            BitSet bits2 = new BitSet();
            
        	BitSet left = new BitSet(); 
        	logCladeCredibility += getLogCladeCredibility(node.getLeft(), left);
        	BitSet right = new BitSet(); 
        	logCladeCredibility += getLogCladeCredibility(node.getRight(), right);
            bits2.or(left);
            bits2.or(right);

//            for (int i=1; i<bits2.length(); i=i+2) {
//                bits2.set(i, false);
//            }

            logCladeCredibility += getLogCladeCredibility(bits2, left, right);

            if (bits != null) {
                bits.or(bits2);
            }
        }

        //System.err.println(node.getNr() + " " + logCladeCredibility);
        return logCladeCredibility;
    }
	
    private double getLogCladeCredibility(BitSet bits, BitSet left, BitSet right) {
    	Clade clade = cladeMap.get(bits);
    	if (clade == null) {
    		return EPSILON;
    	}
        int cladeCount = clade.getCount();
	
    	Map<BitSet, Clade> clades = conditionalCladeMap.get(bits);
        if (clades == null) {
        	return EPSILON;
        }
        
        Clade leftClade = clades.get(left);
        int leftCount = (leftClade != null ? leftClade.getCount() : 0);

        Clade rightClade = clades.get(right);
        int rightCount = (rightClade != null ? rightClade.getCount() : 0);
        
        int subCladCount = Math.max(leftCount, rightCount);
        if (subCladCount == 0) {
        	return EPSILON;
        }
        double logP = Math.log(subCladCount) - Math.log(cladeCount);
        return logP;
    }

	@Override
	public List<String> getArguments() {return null;}

	@Override
	public List<String> getConditions() {return null;}

	@Override
	public void sample(State state, Random random) {}

    public class Clade {
        public Clade(BitSet bits) {
            this.bits = bits;
            count = 0;
            credibility = 0.0;
        }

        public int getCount() {
            return count;
        }

        public void setCount(int count) {
            this.count = count;
        }

        public double getCredibility() {
            return credibility;
        }

        public void setCredibility(double credibility) {
            this.credibility = credibility;
        }

        public List<Object[]> getAttributeValues() {
            return attributeValues;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final Clade clade = (Clade) o;

            return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

        }

        @Override
        public int hashCode() {
            return (bits != null ? bits.hashCode() : 0);
        }

        @Override
        public String toString() {
        	if (tree != null) {
        		List<String> taxa = tree.getTaxonset().asStringList();
        		StringBuilder b = new StringBuilder();
        		b.append("{");
        		for (int i = 0; i < taxa.size(); i++) {
        			if (bits.get(i * 2)) {
        				b.append(taxa.get((int)i)).append(",");
        			}
        		}
        		b.deleteCharAt(b.length() - 1);
        		b.append('}');
        		return b.toString();
        	}
            return "clade " + bits.toString();
        }

        int count;
        double credibility;
        BitSet bits;
        List<Object[]> attributeValues = null;
    }

	public void listConditionalCladeProbabilities() {
		for (BitSet bitset : conditionalCladeMap.keySet()) {
			Map<BitSet, Clade> clades = conditionalCladeMap.get(bitset);
	    	Clade clade = cladeMap.get(bitset);
	        int cladeCount = clade.getCount();

			
			for (BitSet left : clades.keySet()) {
				BitSet right = new BitSet();
				right.or(bitset);
				right.xor(left);
				
		        Clade leftClade = clades.get(left);
		        int leftCount = (leftClade != null ? leftClade.getCount() : 0);

		        Clade rightClade = clades.get(right);
		        int rightCount = (rightClade != null ? rightClade.getCount() : 0);
		        
		        int subCladCount = Math.max(leftCount, rightCount);
		        System.out.print(leftClade + " " + rightClade + " :");
		        if (subCladCount == 0) {
		        	System.out.println(EPSILON);
		        } else {
		        	System.out.println((double)subCladCount/cladeCount);
		        }
			}
		}		
	}

	public void setTree(Tree tree) {
		((Tree) this.tree).assignFrom(tree);
	}
}
