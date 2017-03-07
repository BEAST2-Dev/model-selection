package beast.gss.distributions;

import java.io.IOException;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import beast.app.treeannotator.TreeAnnotator;
import beast.app.treeannotator.TreeAnnotator.TreeSet;
import beast.app.util.TreeFile;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Param;
import beast.core.State;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

@Description("Conditional Clade Distribution for a tree set")
public class CCProbability extends Distribution {
	private TreeFile treeFile;
	private TreeInterface tree;
	private int burninPercentage;
	
	public TreeFile getTreeFile() {return treeFile;}
	public void setTreeFile(TreeFile treeFile) {this.treeFile = treeFile;}
	public TreeInterface getTree() {return tree;}
	public void setTree(TreeInterface tree) {this.tree = tree;}
	
    public int getBurnine() {
		return burninPercentage;
	}
	public void setBurnin(int burninPercentage) {
		this.burninPercentage = burninPercentage;
	}

	protected Map<BitSet, Map<BitSet, Clade>> conditionalCladeMap = new HashMap<>();
    protected Map<BitSet, Clade> cladeMap = new HashMap<>();

    // log probability for a clade that does not exist in the clade system
    final static double EPSILON = -1e100;
	
	public CCProbability(@Param(name="treefile", description="file containing tree set") TreeFile treeFile,
			@Param(name="tree", description="beast tree for which the conditional clade distribution is calculated") TreeInterface tree,
			@Param(name="burnin", description="percentage of the tree set to remove from the beginning") int burninPercentage) {
		this.treeFile = treeFile;
		this.tree = tree;
		this.burninPercentage = burninPercentage;
		if (burninPercentage < 0 || burninPercentage >= 100) {
			throw new IllegalArgumentException("burnin must be a positive number not larger than 100");
		}
		processTreeFile();
	}

	@Override
	public void initAndValidate() {
		super.initAndValidate();
	}
	
	private void processTreeFile() {
		
		try {
			TreeSet trees = new TreeAnnotator().new MemoryFriendlyTreeSet(treeFile.getAbsolutePath(), burninPercentage);
			trees.reset();
			while (trees.hasNext()) {
				Tree tree = trees.next();
				addClades(tree.getRoot());
			}
		} catch (IOException e) {
			e.printStackTrace();
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
	
    private BitSet addClades(Node node) {

        BitSet bits = new BitSet();

        if (node.isLeaf()) {

            int index = node.getNr();
            bits.set(2*index);

        } else {

        	BitSet left = addClades(node.getLeft());
        	BitSet right = addClades(node.getRight());
            bits.or(left);
            bits.or(right);

            for (int i=1; i<bits.length(); i=i+2) {
                bits.set(i, false);
            }
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
		return logP;
	}
	
    public double getLogCladeCredibility(Node node, BitSet bits) {

        double logCladeCredibility = 0.0;

        if (node.isLeaf()) {
            int index = node.getNr();
            bits.set(2*index);
        } else {

            BitSet bits2 = new BitSet();
            
        	BitSet left = new BitSet(); 
        	logCladeCredibility += getLogCladeCredibility(node.getLeft(), left);
        	BitSet right = new BitSet(); 
        	logCladeCredibility += getLogCladeCredibility(node.getRight(), right);
            bits2.or(left);
            bits2.or(right);

            for (int i=1; i<bits2.length(); i=i+2) {
                bits2.set(i, false);
            }

            logCladeCredibility += Math.log(getCladeCredibility(bits2, left, right));

            if (bits != null) {
                bits.or(bits2);
            }
        }

        return logCladeCredibility;
    }
	
    private double getCladeCredibility(BitSet bits, BitSet left, BitSet right) {
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
        
        return Math.log(Math.max(leftCount, rightCount)) - Math.log(cladeCount);
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
            return "clade " + bits.toString();
        }

        int count;
        double credibility;
        BitSet bits;
        List<Object[]> attributeValues = null;
    }
}
