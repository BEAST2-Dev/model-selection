package modelselection.gss;

import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.evolution.tree.Tree;
import modelselection.gss.distribution.GSSTreeDistribution;
import beast.base.core.Input.Validate;

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
