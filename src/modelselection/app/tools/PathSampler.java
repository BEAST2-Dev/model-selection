package modelselection.app.tools;


import beast.app.beauti.Beauti;
import beast.app.beauti.BeautiConfig;
import beast.app.beauti.BeautiDoc;
import beast.app.draw.BEASTObjectDialog;
import beast.app.draw.BEASTObjectPanel;
import beast.app.util.Application;
import beast.app.util.ConsoleApp;
import beast.app.util.XMLFile;
import jam.util.IconUtils;
import modelselection.inference.PathSamplerFromFile;

//command line interface to PathSampler
public class PathSampler {
	
	private static ConsoleApp app;
	
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			// create the class with application that we want to launch
			PathSamplerFromFile sampler = new PathSamplerFromFile();
			
			if (args.length == 0) {
				// try the GUI version
				
				// need to set the ID of the BEAST-object
				sampler.setID("GSS");
				
				// then initialise
				sampler.initAndValidate();
				
				// create BeautiDoc and beauti configuration
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
				
				// suppress a few inputs that we don't want to expose to the user
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".mcmc");
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".value");
				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".hosts");
				
				// check wheter the model1Input is correctly set up
				String fileSep = System.getProperty("file.separator");
				if (!sampler.model1Input.get().exists()) {
					sampler.model1Input.setValue(new XMLFile(Beauti.g_sDir + fileSep + "model.xml"), sampler);
				}
			
				// create panel with entries for the application
				BEASTObjectPanel panel = new BEASTObjectPanel(sampler, sampler.getClass(), doc);
				
				// wrap panel in a dialog
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);

				// show the dialog
				if (dialog.showDialog()) {
					dialog.accept(sampler, doc);
					// create a console to show standard error and standard output
					app = new ConsoleApp("PathSampler", 
							"Path Sampler: " + sampler.model1Input.get().getPath(),
					        IconUtils.getIcon(modelselection.app.tools.PathSampleAnalyser.class, "ps.png"));
					sampler.initAndValidate();
					sampler.run();
				}
				return;
			}

			// continue with the command line version
			main = new Application(sampler);
			main.parseArgs(args, false);
			sampler.initAndValidate();
			sampler.run();
			
		} catch (Exception e) {
			// error handling
			System.out.println(e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
