package modelselection.app.tools;



import beastfx.app.beauti.Beauti;
import beastfx.app.inputeditor.BeautiConfig;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BEASTObjectDialog;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beastfx.app.tools.Application;
import beast.base.core.ProgramStatus;
import beast.pkgmgmt.Utils6;
import beastfx.app.util.XMLFile;
import modelselection.gss.GSSFromFile;

//command line interface to GSS
public class GeneralisedSteppingStone {
	
	public static void main(final String[] args) throws Exception {
		
		// suppress a few inputs that we don't want to expose to the user
		String className = GSSFromFile.class.getName();
		String [] suppressedInputs = new String [] {
				className + ".mcmc",
				className + ".value",
				className + ".hosts"};
		new Application(new GSSFromFile(), suppressedInputs, "Generalised Stepping Stone", args);

//		Application main = null;
//		try {
//			// create the class with application that we want to launch
//			GSSFromFile sampler = new GSSFromFile();
//			
//			if (args.length == 0) {
//				// try the GUI version
//				Utils6.startSplashScreen();
//
//				// need to set the ID of the BEAST-object
//				sampler.setID("Generalised Stepping Stone");
//				
//				// then initialise
//				sampler.initAndValidate();
//				
//				// create BeautiDoc and beauti configuration
//				BeautiDoc doc = new BeautiDoc();
//				doc.beautiConfig = new BeautiConfig();
//				doc.beautiConfig.initAndValidate();
//				
//				// suppress a few inputs that we don't want to expose to the user
//				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".mcmc");
//				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".value");
//				doc.beautiConfig.suppressBEASTObjects.add(sampler.getClass().getName() + ".hosts");
//				
//				// check wheter the model1Input is correctly set up
//				String fileSep = System.getProperty("file.separator");
//				if (!sampler.model1Input.get().exists()) {
//					sampler.model1Input.setValue(new XMLFile(ProgramStatus.g_sDir + fileSep + "model.xml"), sampler);
//				}
//			
//				// create panel with entries for the application
//				BEASTObjectPanel panel = new BEASTObjectPanel(sampler, sampler.getClass(), doc);
//				
//				// wrap panel in a dialog
//				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
//
//				Utils6.endSplashScreen();
//				// show the dialog
//				if (dialog.showDialog()) {
//					dialog.accept(sampler, doc);
//					// create a console to show standard error and standard output
//					app = new ConsoleApp("Generalised Stepping Stone", 
//							"Generalised Stepping Stone: " + sampler.model1Input.get().getPath(),
//					        IconUtils.getIcon(modelselection.app.tools.PathSampleAnalyser.class, "ps.png"));
//					sampler.initAndValidate();
//					sampler.run();
//				}
//				return;
//			}
//
//			// continue with the command line version
//			main = new Application(sampler);
//			main.parseArgs(args, false);
//			sampler.initAndValidate();
//			sampler.run();
//			
//		} catch (Exception e) {
//			// error handling
//			System.out.println(e.getMessage());
//			if (main != null) {
//				System.out.println(main.getUsage());
//			}
//		}
	}

}
