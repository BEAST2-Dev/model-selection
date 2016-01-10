package beast.app.tools;


import java.io.File;

import beast.app.beauti.Beauti;
import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;

// command line interface to PairedPathSampler
public class PairedPathSampler {

	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			beast.inference.PairedPathSampler sampler = new beast.inference.PairedPathSampler();
			
			if (args.length == 0) {
				// try the GUI version
				sampler.setID("PairedPathSampler");
				sampler.initAndValidate();
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".mcmc");
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".value");
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".hosts");
				String fileSep = System.getProperty("file.separator");
				if (!sampler.model1Input.get().exists()) {
					sampler.model1Input.setValue(new File(Beauti.g_sDir + fileSep + "model1.xml"), sampler);
					sampler.model2Input.setValue(new File(Beauti.g_sDir + fileSep + "model2.xml"), sampler);
				}
				
				BEASTObjectPanel panel = new BEASTObjectPanel(sampler, sampler.getClass(), doc);
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
				if (dialog.showDialog()) {
					dialog.accept(sampler, doc);
					sampler.initAndValidate();
					sampler.run();
				}
				System.exit(0);
				return;
			}

			// continue with the command line version
			main = new Application(sampler);
			main.parseArgs(args, false);
			sampler.initAndValidate();
			sampler.run();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
