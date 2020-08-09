package modelselection.gss;

import java.util.List;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.core.StateNodeInitialiser;
import beast.evolution.tree.Tree;
import modelselection.gss.distribution.GSSTreeDistribution;
import beast.core.Input.Validate;

@Description("Initialises tree by randomly picking one from a tree file")
public class TreeFromTreeSetFileInitialiser extends BEASTObject implements StateNodeInitialiser {
	final public Input<GSSTreeDistribution> treeFileInput = new Input<>("treeFile", "file to randomly select starting tree from", Validate.REQUIRED);
	final public Input<Tree> initialInput = new Input<>("initial", "the tree to initialise", Validate.REQUIRED);

	@Override
	public void initAndValidate() {
	}

	@Override
	public void initStateNodes() {
		GSSTreeDistribution treeFile = treeFileInput.get();
		Tree tree = treeFile.getRandomTree();
		initialInput.get().assignFromWithoutID(tree);
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		stateNodes.add(initialInput.get());
	}

}
