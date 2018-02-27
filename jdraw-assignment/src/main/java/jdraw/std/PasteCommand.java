/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import jdraw.framework.DrawCommand;
import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

public class PasteCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private Figure clone;
	private DrawContext context;

	public PasteCommand(Figure clone, DrawContext context) {
		this.clone = clone;
		this.context = context;
	}

	@Override
	public void redo() {
		context.getModel().addFigure(clone);		
		context.getView().addToSelection(clone);
	}

	@Override
	public void undo() {
		context.getModel().removeFigure(clone);		
		context.getView().removeFromSelection(clone);
		ClipBoard.get().get(ClipBoard.get().size()-1).move(-25, -25);	
	}
	
}

