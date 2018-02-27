/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import jdraw.framework.DrawCommand;
import jdraw.framework.DrawModel;
import jdraw.framework.Figure;


public class AddFigureCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private Figure f;

	private DrawModel model;

	public AddFigureCommand(Figure aFigure, DrawModel model) {
		this.f = aFigure;
		this.model = model;
	}


	@Override
	public void redo() {
		model.addFigure(f);
	}

	
	@Override
	public void undo() {
		model.removeFigure(f);
	}

}

