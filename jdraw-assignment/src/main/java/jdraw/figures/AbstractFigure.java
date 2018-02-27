package jdraw.figures;


import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import jdraw.framework.Figure;
import jdraw.framework.FigureEvent;
import jdraw.framework.FigureHandle;
import jdraw.framework.FigureListener;

public abstract class AbstractFigure implements Cloneable, Figure{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	private List<FigureListener> listeners = new CopyOnWriteArrayList <>();
	
	@Override
	public abstract void draw(Graphics g);

	
	@Override
	public abstract void move(int dx, int dy);

	@Override
	public abstract boolean contains(int x, int y);
	
	@Override
	public abstract void setBounds(Point origin, Point corner);

	@Override
	public abstract Rectangle getBounds();

	/**
	 * Returns a list of 8 handles for this Rectangle.
	 * @return all handles that are attached to the targeted figure.
	 * @see jdraw.framework.Figure#getHandles()
	 */	
	@Override
	public abstract List<FigureHandle> getHandles();

	@Override
	public void addFigureListener(FigureListener listener) {
		listeners.add(listener);
	}

	@Override
	public void removeFigureListener(FigureListener listener) {
		listeners.remove(listener);
	}

	protected void propagateFigureEvent() {
		FigureEvent event = new FigureEvent(this);
		for (FigureListener l : listeners) { l.figureChanged(event); }
	}
	
	@Override
	public AbstractFigure clone(){
		try{
			AbstractFigure copy = (AbstractFigure) super.clone();
			copy.listeners = new LinkedList<FigureListener>();
			return copy;
		}catch(CloneNotSupportedException e){
			throw new InternalError("Ciaone");
		}
	}
	
}
