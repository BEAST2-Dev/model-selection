package test.simstudy;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beastfx.app.tools.LogAnalyser;
import beast.base.parser.XMLParserException;
import modelselection.cpo.CPOAnalyser;
import modelselection.cpo.CPOAnalyser.CPOTable;
import modelselection.inference.AICMAnalyser;
import modelselection.inference.AICMAnalyser.AnalysisType;

public class Analyse {

	public static void main(String[] args) throws Exception {
		Analyse a = new Analyse();
		
		for (AnalysisType type : AnalysisType.values()) {
			a.analyse("cc", type);
			a.analyse("ec0.01", type);
			a.analyse("ec0.025", type);
			a.analyse("ec0.05", type);
			a.analyse("ec0.1", type);
		}
		
		//a.cpoAnalysis("cc");
	}

	
	
	private void cpoAnalysis(String dir) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		int n = 0;
		for (int i = 0; i < 100; i++) {
			copyFile(dir +"/cc/run1/cc" + i + ".log", dir +"/cc/cc.log");
			copyFile(dir +"/cc/run1/cc" + i + ".trees", dir +"/cc/cc.trees");
			double score1 = cpo(dir +"/cc/analysis-out" + i + ".xml");

			copyFile(dir +"/ec/run1/ec" + i + ".log", dir +"/ec/ec.log");
			copyFile(dir +"/ec/run1/ec" + i + ".trees", dir +"/ec/ec.trees");
			double score2 = cpo(dir +"/ec/analysis-out" + i + ".xml");
			if (score2 > score1) {
				n++;
			}
		}
		System.out.println(dir + " LPML " + (100-n) + " " + n);
	}
	
    static void copyFile(String source, String target) {
        try {
            Files.copy(new File(source).toPath(), new File(target).toPath());
        } catch (IOException x) {
            System.err.format("Unable to copy: %s: %s%n", source, x);
        }
    }

	
	private double cpo(String file) throws SAXException, IOException, ParserConfigurationException, XMLParserException {
		System.setProperty("java.only", "true");
		CPOAnalyser analyser = new CPOAnalyser();
		analyser.xmlFileInput.setValue(file, analyser);
		
		CPOTable cpoTable = null;
		cpoTable = analyser.getCPOTableFromXML();
		if (cpoTable == null) {
			throw new RuntimeException("Could not get CPOTable from XML");
		}
		
		double [][] patterLogProbs = cpoTable.getPatternLogProbs();
		//int [] weights = cpoTable.getPatternWeights();
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
		double LPML = analyser.calcLPML(order, minLogP, cpoTable);
		return LPML;
	}



	private void analyse(String dir, AICMAnalyser.AnalysisType type) throws IOException {
		AICMAnalyser analyser = new AICMAnalyser();
		analyser.typeInput.setValue(type, analyser);
		analyser.bootstrapLengthInput.setValue(-1, analyser);
		
		int n = 0;
		for (int i = 0; i < 100; i++) {
			double score1 = score(dir +"/cc/run1/cc" + i + ".log", analyser);
			double score2 = score(dir +"/ec/run1/ec" + i + ".log", analyser);
			if (score2 > score1) {
				n++;
			}
		}
		System.out.println(dir + " " + type +  " " + (100-n) + " " + n);
	}

	double score(String file, AICMAnalyser analyser) throws IOException {
		LogAnalyser trace1 = new LogAnalyser(file, 10, true, false);
		Double [] likelihood1 = trace1.getTrace("likelihood");
		double score = analyser.calculateLogMarginalLikelihood(likelihood1);
		return score;
	}


}
