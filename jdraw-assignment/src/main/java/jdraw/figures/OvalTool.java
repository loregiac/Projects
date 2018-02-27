/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.figures;

import java.awt.Point;
import jdraw.framework.DrawContext;
import jdraw.framework.Figure;

/**
 * This tool defines a mode for drawing rectangles.
 *
 * @see jdraw.framework.Figure
 *
 * @author  Christoph Denzler
 */
public class OvalTool extends AbstractDragDrawTool {
	/**
	 * Create a new rectangle tool for the given context.
	 * @param context a context to use this tool in.
	 */
	
	public OvalTool(DrawContext context, String name, String icon) {
		super(context,name,icon);
	}


	@Override
	protected Figure createFigure(Point p) {
		// TODO Auto-generated method stub
		return new Oval(p);
	}

}
