package jdraw.std;

import java.util.Stack;

import jdraw.framework.DrawCommand;

public class CompositeDrawCommand implements DrawCommand {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private Stack<DrawCommand> commands = new Stack<DrawCommand>();
	
	@Override
	public void redo() {
		for(DrawCommand cmd : commands){
			cmd.redo();
		}
	}

	@Override
	public void undo() {
		for(DrawCommand cmd : commands){
			cmd.undo();
		}
	}
	
	public void addCommand(DrawCommand cmd){
		commands.add(cmd);
	}
	
	public Stack<DrawCommand> getScript(){
		return commands;
	}
	
	public void setScript(Stack<DrawCommand> commands){
		this.commands = commands;
	}

}
