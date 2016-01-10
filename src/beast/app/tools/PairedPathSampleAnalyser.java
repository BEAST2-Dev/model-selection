package beast.app.tools;

import java.io.File;

import beast.app.beauti.Beauti;
import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.core.BEASTObject;
import beast.core.Input;
import beast.inference.PathSamplerFromFile;

//command line interface to PathSampler
public class PairedPathSampleAnalyser {
		
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			beast.inference.PairedPathSampleAnalyser analyser = new beast.inference.PairedPathSampleAnalyser();
			
			if (args.length == 0) {
				// try the GUI version
				analyser.setID("PathSamplerAnalyser");
				analyser.initAndValidate();
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
				doc.beautiConfig.suppressBEASTObjects.add(analyser.getClass().getName() + ".mcmc");
				doc.beautiConfig.suppressBEASTObjects.add(analyser.getClass().getName() + ".value");
				doc.beautiConfig.suppressBEASTObjects.add(analyser.getClass().getName() + ".hosts");
				String fileSep = System.getProperty("file.separator");
			
				BEASTObjectPanel panel = new BEASTObjectPanel(analyser, analyser.getClass(), doc);
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
				if (dialog.showDialog()) {
					dialog.accept(analyser, doc);
					analyser.initAndValidate();
					analyser.run();
				}
				System.exit(0);
				return;
			}

			// continue with the command line version
			main = new Application(analyser);
			main.parseArgs(args, false);
			analyser.initAndValidate();
			analyser.run();
		} catch (Exception e) {
			//e.printStackTrace();
			System.out.println(e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
