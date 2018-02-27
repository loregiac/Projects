package jdraw.figures;

import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.MouseEvent;

import jdraw.decorators.DecoratorFigure;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;
import jdraw.framework.FigureHandle;

public class DecoratorHandle implements FigureHandle {

	private FigureHandle handle;
	private DecoratorFigure parent;
	
	public DecoratorHandle(FigureHandle handle, DecoratorFigure parent){
		this.handle = handle;
		this.parent = parent;
	}
	@Override
	public Figure getOwner() {
		return this.parent;
	}

	@Override
	public Point getLocation() {
		return this.handle.getLocation();
	}

	@Override
	public void draw(Graphics g) {
		this.handle.draw(g);
	}

	@Override
	public Cursor getCursor() {
		return this.handle.getCursor();
	}

	@Override
	public boolean contains(int x, int y) {
		return this.handle.contains(x, y);
	}

	@Override
	public void startInteraction(int x, int y, MouseEvent e, DrawView v) {
		this.handle.startInteraction(x, y, e, v);
	}

	@Override
	public void dragInteraction(int x, int y, MouseEvent e, DrawView v) {
		this.handle.dragInteraction(x, y, e, v);
	}

	@Override
	public void stopInteraction(int x, int y, MouseEvent e, DrawView v) {
		this.handle.startInteraction(x, y, e, v);
	}

	@Override
	public HandleState getState() {
		return this.handle.getState();
	}

	@Override
	public void setState(HandleState Nstate) {
		this.handle.setState(Nstate);
	}

}
