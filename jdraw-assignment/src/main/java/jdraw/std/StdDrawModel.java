/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import java.util.LinkedList;


import jdraw.framework.DrawCommandHandler;
import jdraw.framework.DrawModel;
import jdraw.framework.DrawModelEvent;
import jdraw.framework.DrawModelListener;
import jdraw.framework.Figure;
import jdraw.framework.FigureEvent;
import jdraw.framework.FigureListener;

/**
 * Provide a standard behavior for the drawing model. This class initially does not implement the methods
 * in a proper way.
 * It is part of the course assignments to do so.
 * @author TODO add your name here
 *
 */
public class StdDrawModel implements DrawModel, FigureListener {

	private LinkedList<Figure> figures = new LinkedList<Figure>();
	private LinkedList<DrawModelListener> listeners = new LinkedList<DrawModelListener>();
	
	
	@Override
	public void figureChanged(FigureEvent e){
		Figure f = e.getFigure();
		notifyListeners(f,DrawModelEvent.Type.FIGURE_CHANGED);
	}

	public void notifyListeners(Figure f, DrawModelEvent.Type type){
		DrawModelEvent e = new DrawModelEvent(this,f,type);
		for(DrawModelListener l:listeners){
			l.modelChanged(e);
		}
	}
	
	@Override
	public void addFigure(Figure f) {
		// TODO to be implemented
		//System.out.println("StdDrawModel.addFigure has to be implemented");
		//proceed only if the figure is not already contained in the list
		if(!figures.contains(f)){
			// add FigureListener
			figures.add(f);
			notifyListeners(f, DrawModelEvent.Type.FIGURE_ADDED);
			f.addFigureListener(this);
		}
	}
	
	 
	
	@Override
	public Iterable<Figure> getFigures() {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.getFigures has to be implemented");
		//return new LinkedList<Figure>(); // Only guarantees, that the application starts -- has to be replaced !!!
		return figures; 
		
	}

	@Override
	public void removeFigure(Figure f) {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.removeFigure has to be implemented");
		
		if(figures.contains(f)){
			//handler.addCommand(new RemoveFigureCommand(this,f));
			figures.remove(f);
			notifyListeners(f, DrawModelEvent.Type.FIGURE_REMOVED);
			f.removeFigureListener(this);
		}
	}

	@Override
	public void addModelChangeListener(DrawModelListener listener) {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.addModelChangeListener has to be implemented");
		listeners.add(listener);
	}

	@Override
	public void removeModelChangeListener(DrawModelListener listener) {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.removeModelChangeListener has to be implemented");
		listeners.remove(listener);
	}

	/** The draw command handler. Initialized here with a dummy implementation. */
	// TODO initialize with your implementation of the undo/redo-assignment.
	//private DrawCommandHandler handler = new EmptyDrawCommandHandler();
	private DrawCommandHandler handler = new ListDrawCommandHandler();


	/**
	 * Retrieve the draw command handler in use.
	 * @return the draw command handler.
	 */
	@Override
	public DrawCommandHandler getDrawCommandHandler() {
		return handler;
	}

	@Override
	public void setFigureIndex(Figure f, int index) {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.setFigureIndex has to be implemented");
		if(!figures.contains(f)){
			throw new IllegalArgumentException("Input figure is not contained in the model");
		}
		if(index > figures.size()-1){
			throw new IndexOutOfBoundsException("Input index is out of range");
		}
		int old = figures.indexOf(f);
		figures.add(index, figures.remove(old));
		notifyListeners(f, DrawModelEvent.Type.DRAWING_CHANGED);
		
	}

	@Override
	public void removeAllFigures() {
		// TODO to be implemented  
		//System.out.println("StdDrawModel.removeAllFigures has to be implemented");
		for(Figure f : figures){
			f.removeFigureListener(this);
		}
		figures = new LinkedList<Figure>();
		notifyListeners(null, DrawModelEvent.Type.DRAWING_CLEARED);

	}

}


