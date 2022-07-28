package modelselection.app.tools;

import beastfx.app.inputeditor.BeautiConfig;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BEASTObjectDialog;
import beastfx.app.inputeditor.BEASTObjectPanel;
import beastfx.app.tools.Application;
import beastfx.app.util.ConsoleApp;
import jam.util.IconUtils;

//command line interface to PathSampleAnalyser
public class PathSampleAnalyser {
		
	public static void main(final String[] args) throws Exception {
		Application main = null;
		try {
			// create the runnable class with application that we want to launch
			modelselection.inference.PathSampleAnalyser analyser = new modelselection.inference.PathSampleAnalyser();
			
			if (args.length == 0) {
				// try the GUI version

				// need to set the ID of the BEAST-object
				analyser.setID("PathSampleAnalyser");
				
				// then initialise
				analyser.initAndValidate();
				
				// create BeautiDoc and beauti configuration
				BeautiDoc doc = new BeautiDoc();
				doc.beautiConfig = new BeautiConfig();
				doc.beautiConfig.initAndValidate();
			
				// create panel with entries for the application
				BEASTObjectPanel panel = new BEASTObjectPanel(analyser, analyser.getClass(), doc);
				
				// wrap panel in a dialog
				BEASTObjectDialog dialog = new BEASTObjectDialog(panel, null);
				if (dialog.showDialog()) {
					dialog.accept(analyser, doc);
					analyser.initAndValidate();

					// create a console to show standard error and standard output
					analyser.consoleApp = new ConsoleApp("PathSampleAnalyser", // name 
							"Path Sample Analyser: " + analyser.rootDirInput.get(), // console title
					        IconUtils.getIcon(modelselection.app.tools.PathSampleAnalyser.class, "ps.png")
							);

					analyser.run();
				}
				return;
			}

			// continue with the command line version
			main = new Application(analyser);
			main.parseArgs(args, false);
			analyser.initAndValidate();
			analyser.run();
		} catch (Exception e) {
			System.out.println(e.getMessage());
			if (main != null) {
				System.out.println(main.getUsage());
			}
		}
	}

}
