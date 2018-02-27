/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import java.util.Stack;

import jdraw.framework.DrawCommand;
import jdraw.framework.DrawCommandHandler;

public class ListDrawCommandHandler implements DrawCommandHandler {

	private Stack<DrawCommand> stack = new Stack<DrawCommand>();
	private int present;
	private boolean scriptOpen = false;
	private CompositeDrawCommand script;
	
	@Override
	public void addCommand(DrawCommand cmd) {
		if(!scriptOpen){
			stack.add(cmd);
			present = stack.size();
		}else{
			script.addCommand(cmd);
		}
	}
	
	@Override
	public void undo() { 
		stack.get(present-1).undo();
		present = present-1;
	}

	@Override
	public void redo() { 
		stack.get(present).redo();
		present = present +1;
	}

	@Override
	public boolean undoPossible() { 
		return present > 0;
	}

	@Override
	public boolean redoPossible() { 
		return present < stack.size(); 
	}

	@Override
	public void beginScript() { 
		scriptOpen = true;
		script = new CompositeDrawCommand();
	}

	@Override
	public void endScript() { 
		scriptOpen = false;
		stack.add(script);
		present = stack.size();
	}

	@Override
	public void clearHistory() { 
		stack.clear(); 
	}
}
