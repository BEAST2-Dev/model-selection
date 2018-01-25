package test.beast.gss.distributions;

import java.io.FileNotFoundException;
import java.io.PrintStream;

import org.junit.Test;

import beast.app.util.TreeFile;
import beast.gss.distributions.GSSTreeDistribution;
import beast.util.TreeParser;
import junit.framework.TestCase;

public class GSSTreeDistributionTest extends TestCase {
	
	@Test
	public void testGSSTreeDistribution() throws FileNotFoundException {
		TreeParser tree = new TreeParser("((A:0.1,B:0.1):0.1,C:0.2)");
		TreeFile file = new TreeFile("/tmp/gsstest.trees");
        PrintStream out = new PrintStream(file.getAbsolutePath());
        tree.init(out);
        out.println();
        tree.log(0, out);
        out.println();
        tree.getRoot().setHeight(0.3);
        tree.getRoot().getLeft().setHeight(0.2);
        tree.log(0, out);
        out.println();        
        tree.close(out);
        out.close();

		GSSTreeDistribution distr  = new GSSTreeDistribution();
		distr.initByName("treefile", file, "tree", tree, "burnin", 0, "useGammaForBranchLengths", GSSTreeDistribution.BranchLengthDistribution.useIntervals);
		double logP = distr.calculateLogP();
		assertEquals(4.161503004629573, logP, 1e-10);
				
	}

	@Test
	public void testGSSTreeDistribution2() {
		GSSTreeDistribution distr  = new GSSTreeDistribution();
		distr.initByName("treefile", "/tmp/one.trees", "tree", null, "burnin", 0, "useGammaForBranchLengths", GSSTreeDistribution.BranchLengthDistribution.none);
		distr.listConditionalCladeProbabilities();
		
	}
	
}
