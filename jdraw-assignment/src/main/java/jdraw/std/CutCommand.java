/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import jdraw.framework.DrawCommand;
import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

public class CutCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	private Figure figure;
	private DrawContext context;
	private Figure clone;

	public CutCommand(Figure clone,Figure figure, DrawContext context) {
		this.figure = figure;
		this.context = context;
		this.clone = clone;
	}

	@Override
	public void redo() {
		ClipBoard.add(clone);
		context.getView().removeFromSelection(figure);
		context.getModel().removeFigure(figure);
	}

	@Override
	public void undo() {
		ClipBoard.remove(clone);
		context.getModel().addFigure(figure);
		context.getView().addToSelection(figure);
	}

}

