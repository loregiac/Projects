package jdraw.decorators;

import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import jdraw.figures.AbstractFigure;
import jdraw.figures.HandleState;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;
import jdraw.framework.FigureEvent;
import jdraw.framework.FigureHandle;
import jdraw.framework.FigureListener;

public class DecoratorFigure extends AbstractFigure implements FigureListener{
	int dx, dy;
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private Figure inner;
	
	public DecoratorFigure(Figure figure){
		inner = figure;
		inner.addFigureListener(this);
	}
	
	@Override
	public void swapHorizontal() {
		inner.swapHorizontal();
	}

	@Override
	public void swapVertical() {
		inner.swapVertical();
	}

	@Override
	public void draw(Graphics g) {
		inner.draw(g);
	}

	@Override
	public void move(int dx, int dy) {
		this.dx = dx;
		this.dy = dy;
		inner.move(dx, dy);
	}

	@Override
	public boolean contains(int x, int y) {
		return inner.contains(x, y);
	}

	@Override
	public void setBounds(Point origin, Point corner) {
		inner.setBounds(origin, corner);
	}

	@Override
	public Rectangle getBounds() {
		return inner.getBounds();
	}

	@Override
	public List<FigureHandle> getHandles() {
		List<FigureHandle> handles = new ArrayList<FigureHandle>(inner.getHandles().size());
		for(FigureHandle h : inner.getHandles()){
			handles.add(new DecoratorHandle(h));
		}
		return Collections.unmodifiableList(handles);
	}

	public Figure getInner(){
		return inner;
	}

	@Override
	public void addFigureListener(FigureListener listener) {
		inner.addFigureListener(listener);
	}

	@Override
	public void removeFigureListener(FigureListener listener) {
		inner.removeFigureListener(listener);
	}

	@Override
	public DecoratorFigure clone() {
		DecoratorFigure copy = (DecoratorFigure) super.clone();
		copy.inner = inner.clone();
		copy.inner.addFigureListener(copy);
		return copy;
	}


	@Override
	public void figureChanged(FigureEvent e) {
		propagateFigureEvent();
	}
	
	public class DecoratorHandle implements FigureHandle {

		private FigureHandle handle;
		
		public DecoratorHandle(FigureHandle handle){
			this.handle = handle;
		}
		@Override
		public Figure getOwner() {
			return DecoratorFigure.this;
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

	@Override
	public int getDx() {
		// TODO Auto-generated method stub
		return dx;
	}

	@Override
	public int getDy() {
		// TODO Auto-generated method stub
		return dy;
	}

}
