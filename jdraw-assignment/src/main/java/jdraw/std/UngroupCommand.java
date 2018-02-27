/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import jdraw.figures.GroupFigure;
import jdraw.framework.DrawCommand;
import jdraw.framework.DrawContext;
import jdraw.framework.DrawModel;
import jdraw.framework.Figure;


public class UngroupCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private List<Figure> groupedFigures = new CopyOnWriteArrayList<>();
	private DrawModel model;
	private GroupFigure group;
	private DrawContext context;

	public UngroupCommand(List<Figure> groupedFigures, GroupFigure group, DrawModel model,DrawContext context) {
		this.groupedFigures = groupedFigures;
		this.group = group;
		this.model = model;
		this.context = context;
	}

	@Override
	public void redo() {
		model.removeFigure(group);
		for(Figure f : groupedFigures){
			model.addFigure(f);
			context.getView().addToSelection(f);
		}
	}

	@Override
	public void undo() {
		for(Figure f : groupedFigures){
			model.removeFigure(f);
		}
		model.addFigure(group);
		context.getView().addToSelection(group);
	}
}