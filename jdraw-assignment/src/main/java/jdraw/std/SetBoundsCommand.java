/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import java.awt.Point;
import java.awt.Rectangle;

import jdraw.framework.DrawCommand;
import jdraw.framework.Figure;

public class SetBoundsCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private Figure f;
	
	private Rectangle oldBounds;
	
	private Rectangle newBounds;

	public SetBoundsCommand(Figure aFigure, Rectangle oldBounds, Rectangle newBounds) {
		this.f = aFigure;
		this.oldBounds = oldBounds;
		this.newBounds = newBounds;
	}


	@Override
	public void redo() {
		Point origin = new Point(newBounds.x,newBounds.y);
		Point corner = new Point(newBounds.x+newBounds.width,newBounds.y+newBounds.height);
		f.setBounds(origin, corner);
	}

	@Override
	public void undo() {
		Point origin = new Point(oldBounds.x,oldBounds.y);
		Point corner = new Point(oldBounds.x+oldBounds.width,oldBounds.y+oldBounds.height);
		f.setBounds(origin, corner);
	}

}

